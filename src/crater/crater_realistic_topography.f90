!***** crater/crater_realistic_topography
! Name
!   crater_realistic_topography 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_crater
!   
!   call crater_realistic_topography(user,surf,crater,domain,deltaMtot)
!
! DESCRIPTION
!   Creates realistic topography for each fresh crater using Perlin noise. 
!    
! 
! ARGUMENTS
!   Input
!   * user    -- User input parameters
!   * surf    -- Surface grid
!   * crater  -- Crater dimension container
!   * domain  -- Simulation domain variable container
!   * deltaMtot -- Total displaced mass used to calculate mass-conserving ejecta blanket
!
!   Output
!   * surf   -- Outputs the new ejecta blanket onto the grid
!   * crater -- May affects the value of the maximum affected distance
! 
! Notes
!
!**********************************************************************************************************************************
subroutine crater_realistic_topography(user,surf,crater,domain,deltaMtot)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: lradsq,xbar, ybar 
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   real(DP) :: r,xp,yp,deltaMi,rimheight,rad
   logical :: lastloop
   real(DP), dimension(2) :: rn

   interface
      subroutine complex_floor(user,surf,crater,rn,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_floor

      subroutine complex_peak(user,surf,crater,rn,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_peak

      subroutine complex_wall_texture(user,surf,crater,rn,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(inout) :: deltaMtot

      end subroutine complex_wall_texture

      subroutine complex_terrace(user,surf,crater,rn,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_terrace

   end interface

   ! Executable code
   call random_number(rn)
   
   call complex_terrace(user,surf,crater,rn,deltaMtot)
   call complex_wall_texture(user,surf,crater,rn,deltaMtot)
   call crater_slope_collapse(user,surf,crater,domain,(0.35_DP * user%pix)**2,deltaMtot)
   call complex_floor(user,surf,crater,rn,deltaMtot)
   call complex_peak(user,surf,crater,rn,deltaMtot)



   return
end subroutine crater_realistic_topography

subroutine complex_peak(user,surf,crater,rn,deltaMtot)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,xbar,ybar,areafrac,r
   integer(I4B) :: i,j,inc,xpi,ypi

   ! Complex crater peak
   real(DP), parameter :: rad = 0.20_DP      ! Radial size of central peak relative to frad
   integer(I4B), parameter :: num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 1000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 16.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 0.80_DP  ! The relative size scaling at each octave level

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise


   inc = max(min(nint(rad * crater%frad / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         areafrac = util_area_intersection(rad * crater%frad,xbar,ybar,user%pix)

         r = sqrt(xbar**2 + ybar**2) / crater%frad

         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
            znoise = crater%peakheight * (1.0_DP - r / rad) * pers ** (octave - 1)
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2)) * znoise
         end do 
         newdem = newdem + max(noise,0.0_DP) * areafrac 

         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem
      end do
   end do
   return

end subroutine complex_peak

subroutine complex_floor(user,surf,crater,rn,deltaMtot)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,xbar,ybar,areafrac
   integer(I4B) :: i,j,inc,xpi,ypi

   ! Complex crater floor parameters
   integer(I4B), parameter :: num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 10000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 4.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: noise_height = 1e-3_DP ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 0.80_DP  ! The relative size scaling at each octave level


   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise

   inc = max(min(nint(0.5_DP * crater%floordiam / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         areafrac = util_area_intersection(0.5_DP * crater%floordiam,xbar,ybar,user%pix)

         ! Make the floor topographically "noisy"
         noise = 0.0_DP
         do octave = 1, num_octaves
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
            znoise = noise_height * pers ** (octave - 1) * crater%fcrat
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2)) * znoise
         end do
         newdem = newdem + max(noise,0.0_DP) * areafrac

         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem
      end do
   end do

   return
end subroutine complex_floor


subroutine complex_wall_texture(user,surf,crater,rn,deltaMtot)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,xbar,ybar,areafrac,flr,r
   integer(I4B) :: i,j,inc,xpi,ypi

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise,dnoise

   ! Complex crater all texture parameters
   real(DP), parameter :: rad = 0.20_DP      ! Radial size of central peak relative to frad
   integer(I4B), parameter :: num_octaves  = 10   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 4000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 4.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: noise_height = 5.0e-3_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 1.20_DP  ! The relative size scaling at each octave level


   inc = max(min(nint(2.0_DP * crater%frad / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   flr = crater%floordiam / crater%fcrat

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         r = sqrt(xbar**2 + ybar**2) / crater%frad

         areafrac = 1.0 - util_area_intersection(0.3_DP * crater%floordiam,xbar,ybar,user%pix)
         areafrac = areafrac * util_area_intersection(2.0_DP * crater%frad,xbar,ybar,user%pix)
         areafrac = areafrac * min(min(r**(-8),1.0_DP),(r / flr)**12,1.0_DP)

         ! Add some roughness to the walls
         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
            znoise = noise_height  * (pers * min(flr / r, r / flr)) ** (octave - 1) * crater%fcrat
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2))* znoise
         end do
         newdem = max(newdem + noise * areafrac,crater%melev - crater%floordepth)

         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem
      end do
   end do

   return
end subroutine complex_wall_texture


subroutine complex_terrace(user,surf,crater,rn,deltaMtot)
   !Makes terraced walls by applying noisy topographic diffusion in discrete radial zones along the wall
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,rad,xbar,ybar,r,flr,hprof,tprof
   integer(I4B) :: xpi,ypi,i,j,inc,maxhits

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   real(DP) :: noise,dnoise

   ! Complex crater wall
   integer(I4B)            :: nscallops             ! Approximate number of scallop features on edge of crater
   integer(I4B),parameter  :: offset = 3000         ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP),parameter      :: noise_height = 0.15_DP ! Vertical height of noise features as a function of crater radius at the first octave
                                                    ! Higher values make the scallop features broader
   real(DP),parameter      :: noise_slope = 0.6_DP  ! The slope of the power law function that defines the noise shape for scallop features
                                                    ! Higher values makes the inner sides of the scallop features more rounded, lower values 
                                                    ! make them have sharper points at the inflections
   integer(I4B),parameter        :: nterraces = 6   ! Number of terraces
   real(DP),parameter            :: terracefac = 2.0_DP ! Power law scaling for relative terrace sizes
   integer(I4B)                  :: terrace         ! The current terrace
   real(DP)                      :: router          ! The radius of the terrace outer edge
   real(DP)                      :: rinner          ! The radius of the terrace outer edge
   real(DP)                      :: upshift,dfloor
   

   ! Internal variables
   real(DP),dimension(:,:),allocatable :: kdiff,cumulative_elchange,areafrac
   integer(I4B),dimension(:,:,:),allocatable :: indarray


   nscallops = 16
   rad = 2.0_DP * crater%frad
   upshift = 0._DP
 
   inc = max(min(nint(rad / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   flr = crater%floordiam / crater%fcrat
   do terrace = 1, nterraces 
      router = flr + (1._DP - flr) * (terrace / real(nterraces,kind=DP))**(terracefac) ! The radius of the outer edge of the terrace
      rinner = flr + (1._DP - flr) * ((terrace - 1) / real(nterraces,kind=DP))**(terracefac) ! The radius of the inner edge of the terrace

      do j = -inc,inc
         do i = -inc,inc
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            newdem = surf(xpi,ypi)%dem

            xbar = xpi * user%pix - crater%xl 
            ybar = ypi * user%pix - crater%yl

            r = sqrt(xbar**2 + ybar**2) / crater%frad

            znoise = noise_height**(1._DP / (2 * noise_slope)) 
            xynoise = (nscallops / PI) / crater%fcrat
            dnoise = util_perlin_noise(xynoise * xbar + terrace * offset * rn(1), &
                                       xynoise * ybar + terrace * offset * rn(2)) * znoise
            noise = (dnoise**2)**noise_slope

            hprof = (r / router)**(-1)
            tprof = crater_profile(user,crater,rinner) * (r/rinner)**2

            if (1.0_DP - noise < hprof) then
               newdem = min(tprof,newdem) 
               elchange  = newdem - surf(xpi,ypi)%dem
               deltaMtot = deltaMtot + elchange
               surf(xpi,ypi)%dem = newdem
               ! Save the minimum elevation below the floor
               dfloor = crater%melev - crater%floordepth - newdem
               if (dfloor > 0.0_DP) upshift = max(upshift,dfloor) 
            end if

         end do
      end do

   end do

   ! Now make the scalloped rim
   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         r = sqrt(xbar**2 + ybar**2) / crater%frad

         ! Make scalloped rim
         znoise = noise_height**(1._DP / (2 * noise_slope))
         xynoise = (nscallops / PI) / crater%fcrat
         dnoise = util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                    xynoise * ybar + offset * rn(2)) * znoise
         noise = (dnoise**2)**noise_slope

         hprof = r**(-1)
         tprof = crater%melev 

         if (1.0_DP - noise < hprof) then
            newdem = min(tprof,newdem) 
            elchange  = newdem - surf(xpi,ypi)%dem
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         end if

      end do
   end do

   ! Compensate for any slumping and lift everything up
   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         r = sqrt(xbar**2 + ybar**2) / crater%frad
         
         if (r > flr) then
            newdem = newdem + upshift * max(min(3.0_DP - 2 * r,1.0_DP),0.0_DP)
            elchange  = newdem - surf(xpi,ypi)%dem
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         end if
      end do 
   end do

   return
end subroutine complex_terrace


