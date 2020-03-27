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

      subroutine complex_rim(user,surf,crater,rn,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_rim
   end interface

   ! Executable code
   call random_number(rn)
   
   call complex_rim(user,surf,crater,rn,deltaMtot)
   call complex_terrace(user,surf,crater,rn,deltaMtot)
   call complex_wall_texture(user,surf,crater,rn,deltaMtot)
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
   real(DP), parameter :: xy_noise_fac = 32.0_DP  ! Spatial "size" of noise features at the first octave
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
   integer(I4B), parameter :: num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 4000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 32.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: noise_height = 1.0e-3_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 0.80_DP  ! The relative size scaling at each octave level


   inc = max(min(nint(1.5_DP * crater%frad / user%pix),PBCLIM*user%gridsize),1) + 1
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
         areafrac = 1.0 - util_area_intersection(0.4_DP * crater%floordiam,xbar,ybar,user%pix)
         areafrac = areafrac * util_area_intersection(1.5_DP * crater%frad,xbar,ybar,user%pix)
         areafrac = areafrac / max(r**2,0.1_DP)

         ! Add some roughness to the walls
         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
            znoise = noise_height  * pers ** (octave - 1) * crater%fcrat
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2))* znoise
         end do
         newdem = newdem + noise * areafrac

         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem
      end do
   end do

   return
end subroutine complex_wall_texture

subroutine complex_rim(user,surf,crater,rn,deltaMtot)
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
   real(DP) :: newdem,elchange,rad,xbar,ybar
   integer(I4B) :: xpi,ypi,i,j,inc,maxhits

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   real(DP) :: noise,dnoise

   ! Complex crater wall
   integer(I4B)        :: num_octaves  ! Number of Perlin noise octaves
   integer(I4B)        :: offset       ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP)            :: xy_noise_fac ! Spatial "size" of noise features at the first octave
   real(DP)            :: freq         ! Spatial size scale factor multiplier at each octave level
   real(DP)            :: pers         ! The relative size scaling at each octave level
   real(DP)            :: noise_height ! Vertical height of noise features as a function of crater radius at the first octave

   ! determine area to effect
   ! First make the interior of the crater
   rad = 1.25_DP * crater%frad
 
   ! Create diffusion noise for terraces
   num_octaves = 4
   offset = 3000
   xy_noise_fac = 5.00_DP
   noise_height = 1.00_DP

   inc = max(min(nint(rad / user%pix),PBCLIM*user%gridsize),1) + 1
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

         ! Make scalloped rim
         znoise = noise_height * crater%fcrat
         xynoise = xy_noise_fac / crater%fcrat
         dnoise = util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                    xynoise * ybar + offset * rn(2)) * znoise
         noise = sqrt(sqrt(dnoise**2))
         !newdem = min(crater%melev + crater%ejrim - noise , newdem) 
         if (crater%melev + crater%ejrim - noise < newdem) newdem = newdem - 1._DP * crater%ejrim

         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem
      end do
   end do

   return

end subroutine complex_rim


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
   real(DP) :: newdem,elchange,rad,xbar,ybar,r,flr
   integer(I4B) :: xpi,ypi,i,j,inc,maxhits

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise,dnoise

   ! Complex crater wall
   integer(I4B)        :: num_octaves  ! Number of Perlin noise octaves
   integer(I4B)        :: offset       ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP)            :: xy_noise_fac ! Spatial "size" of noise features at the first octave
   real(DP)            :: freq         ! Spatial size scale factor multiplier at each octave level
   real(DP)            :: pers         ! The relative size scaling at each octave level
   real(DP)            :: noise_height ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP)            :: scal_height ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP)            :: diff_height  ! Magnitude of diffusion noise 
   integer(I4B)        :: terracefac   ! Spatial size factor for terraces relative to crater radius
   integer(I4B)        :: terrace_num  ! The size of each terrace will be crater radius / terrace_num

   ! Internal variables
   real(DP),dimension(:,:),allocatable :: kdiff,cumulative_elchange,areafrac
   integer(I4B),dimension(:,:,:),allocatable :: indarray


   ! determine area to effect
   ! First make the interior of the crater
   rad = 1.25_DP * crater%frad
 

   ! Create diffusion noise for terraces
   num_octaves = 4
   xy_noise_fac = 4.00_DP
   diff_height = 8e-5_DP * (crater%fcrat)**2
   noise_height = 0.02_DP
   freq = 2.0_DP
   pers = 0.80_DP 
   terracefac = 8
   flr = crater%floordiam / crater%fcrat

   inc = max(min(nint(rad / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   allocate(kdiff(-inc:inc,-inc:inc))
   allocate(areafrac(-inc:inc,-inc:inc))
   allocate(indarray(2,-inc:inc,-inc:inc))

   cumulative_elchange = 0.0_DP
   kdiff = 0.0_DP
   indarray = inc - 1
   areafrac = 0.0_DP

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi

         ! Cut a hole out from the floor and also outside the rim of the crater
         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         r = sqrt(xbar**2 + ybar**2) / crater%frad

         areafrac(i,j) =  (1.0_DP - util_area_intersection(0.95_DP * crater%floordiam / 2,xbar,ybar,user%pix)) 
         areafrac(i,j) = areafrac(i,j) * util_area_intersection(2.0_DP*crater%frad,xbar,ybar,user%pix)

         areafrac(i,j) = areafrac(i,j) / max(r**4,1._DP)

         ! Make the terraces
         ! This shifts the Perlin noise pattern at discrete radial distances, simulating terraces
         terrace_num = int(terracefac * ((r - 1._DP) / (flr - 1._DP))**(1.5_DP)) 
         noise = 0.0_DP
         if ((r <= 1.0_DP).and.(r >= flr)) then
            offset = 2000
            xynoise = xy_noise_fac / crater%fcrat
            znoise  = noise_height * crater%fcrat
            noise = noise + util_perlin_noise(xynoise * xbar + terrace_num * offset * rn(1), &
                                              xynoise * ybar + terrace_num * offset * rn(2))* znoise
            newdem = max(newdem + noise,crater%melev - crater%floordepth)
         end if
         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem

         ! Now smooth out the sharp edges of the terraces with noisy diffusion
         noise = 0.0_DP
         if (r <= 1.0_DP) then
            do octave = 1, num_octaves
               xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
               znoise  = diff_height * pers ** (octave - 1) 
               noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                                 xynoise * ybar + offset * rn(2)) * znoise
            end do
         end if
         kdiff(i,j) = noise 

      end do
   end do
   kdiff(:,:) = kdiff(:,:) - minval(kdiff(:,:))

   ! Cut out the holes
   kdiff(:,:) = kdiff(:,:) * areafrac(:,:) !* 0.0_DP

   maxhits = 1
   call util_diffusion_solver(user,surf,2 * inc + 1,indarray,kdiff,cumulative_elchange,maxhits)

   do j = -inc,inc
      do i = -inc,inc
         xpi = indarray(1,i,j)
         ypi = indarray(2,i,j)
         surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cumulative_elchange(i,j)
         surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cumulative_elchange(i,j),0.0_DP)
      end do
   end do

   deallocate(cumulative_elchange,kdiff,indarray,areafrac)

   return
end subroutine complex_terrace


