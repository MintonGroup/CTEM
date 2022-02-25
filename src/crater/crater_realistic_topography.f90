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
subroutine crater_realistic_topography(user,surf,crater,domain,ejecta_dem)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(inout) :: ejecta_dem

   ! Internal variables
   real(DP) :: deltaMtot
   integer(I4B) :: inc
   real(DP),parameter :: complex_collapse_slope = 0.33_DP  ! Complex craters and basins undergo an extra strong slope collapse

   interface
      subroutine complex_floor(user,surf,crater,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_floor

      subroutine complex_peak(user,surf,crater,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_peak

      subroutine complex_wall_texture(user,surf,crater,domain,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(inout) :: deltaMtot

      end subroutine complex_wall_texture

      subroutine complex_terrace(user,surf,crater,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),intent(inout) :: deltaMtot
      end subroutine complex_terrace

      subroutine ejecta_texture(user,surf,crater,deltaMtot,inc,ejecta_dem)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),intent(in) :: deltaMtot
      integer(I4B),intent(in) :: inc
      real(DP),dimension(-inc:inc,-inc:inc),intent(inout) :: ejecta_dem
      end subroutine ejecta_texture

      subroutine crater_realistic_slope_texture(user,critical_value,inc,critarray)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      real(DP),intent(in) :: critical_value
      integer(I4B),intent(in) :: inc
      real(DP),dimension(-inc:inc,-inc:inc),intent(out) :: critarray
      end subroutine crater_realistic_slope_texture


   end interface

   deltaMtot = 0.0_DP
   select case(crater%morphtype)
   case("COMPLEX","PEAKRING","MULTIRING")
      call complex_terrace(user,surf,crater,deltaMtot)
      call complex_wall_texture(user,surf,crater,domain,deltaMtot)
      call complex_floor(user,surf,crater,deltaMtot)
      call complex_peak(user,surf,crater,deltaMtot)
   end select

   ! Retrieve the size of the ejecta dem and correct for indexing
   inc = (size(ejecta_dem,1) - 1) / 2
   call ejecta_texture(user,surf,crater,deltaMtot,inc,ejecta_dem)

   ! Do a final pass of the slope collapse with a shallower slope than normal to smooth out all of the sharp edges

   if (user%docollapse) then
      select case(crater%morphtype)
      case("COMPLEX","PEAKRING","MULTIRING")
         call crater_slope_collapse(user,surf,crater,domain,(complex_collapse_slope * user%pix)**2,deltaMtot)
      end select
   end if

   return
end subroutine crater_realistic_topography

subroutine complex_peak(user,surf,crater,deltaMtot)
   ! Makes the central peak or peak ring, depending on the crater size
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,xbar,ybar,areafrac,r
   integer(I4B) :: i,j,inc,xpi,ypi
   real(DP), dimension(2) :: rn

   ! Complex crater peak
   integer(I4B), parameter :: num_octaves  = 10   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 1000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 14.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 0.50_DP  ! The relative size scaling at each octave level

   real(DP) :: FWHM,a,b,c,z ! Gaussian parameters

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise


   ! Executable code
   call random_number(rn)

   !Preliminary empirical fits based on a handful of complex craters
   !More work needs to be done to make these models more robust
   FWHM = 0.1_DP !Copernicus
   !FWHM = 0.3_DP !Lansberg 
   a = crater%peakheight
   b = 0.003_DP * ((1e-3_DP * crater%fcrat)**(1.75_DP)) / (1e-3_DP * crater%fcrat) ! Make peak rings for sufficiently large craters
   b = min(b, 0.5_DP)
   c = FWHM / (2 * sqrt(2 * log(2._DP)))
   !*********************

   inc = max(min(nint(0.5_DP * crater%floordiam / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         areafrac = util_area_intersection(0.5_DP * crater%floordiam,xbar,ybar,user%pix)

         r = sqrt(xbar**2 + ybar**2) / crater%frad
         ! Model central peak / peak ring as Gaussian function
         z = a * exp(-(r - b)**2 / (2 * c**2))

         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
            znoise = z * pers ** (octave - 1)
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

subroutine complex_floor(user,surf,crater,deltaMtot)
   ! Makes the flat floor with the appearance of blocks and smooth melts
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,xbar,ybar,areafrac
   integer(I4B) :: i,j,inc,xpi,ypi
   real(DP), dimension(2) :: rn

   ! Complex crater floor parameters
   integer(I4B), parameter :: num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 10000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 4.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: noise_height = 2e-3_DP ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 0.80_DP  ! The relative size scaling at each octave level

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise

   !Executable code
   call random_number(rn)

   inc = max(min(nint(0.5_DP * crater%floordiam / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

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


subroutine complex_wall_texture(user,surf,crater,domain,deltaMtot)
   ! Gives the terraced slopes some rough texture to increase realism
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
   real(DP) :: newdem,elchange,xbar,ybar,areafrac,flr,r
   integer(I4B) :: i,j,inc,xpi,ypi
   real(DP), dimension(2) :: rn

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise,dnoise

   ! Complex crater all texture parameters
   real(DP), parameter :: rad = 0.20_DP      ! Radial size of central peak relative to frad
   integer(I4B), parameter :: num_octaves  = 10   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 4000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 3.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: noise_height = 3.0e-3_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 1.20_DP  ! The relative size scaling at each octave level
   real(DP), parameter :: outer_wall_size = 1.2_DP

   !Executable code
   call random_number(rn)

   inc = max(min(nint(outer_wall_size * crater%frad / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   flr = crater%floordiam / crater%fcrat

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         r = sqrt(xbar**2 + ybar**2) / crater%frad

         areafrac = 1.0 - util_area_intersection(0.3_DP * crater%floordiam,xbar,ybar,user%pix)
         areafrac = areafrac * util_area_intersection(outer_wall_size * crater%frad,xbar,ybar,user%pix)
         areafrac = areafrac * min((r / flr)**12,1.0_DP) ! Smooth out interface between wall and floor
         areafrac = areafrac * max(min(outer_wall_size- r,1.0_DP),0.0_DP) ! Smooth out region outside of the rim

         ! Add some roughness to the walls
         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat
            znoise = noise_height  * (pers * min(flr / r, r / flr)) ** (octave - 1) * crater%fcrat
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2))* znoise
         end do
         newdem = newdem + noise * areafrac
         if (r < flr) newdem = max(newdem,crater%melev - crater%floordepth)

         elchange  = newdem - surf(xpi,ypi)%dem
         deltaMtot = deltaMtot + elchange
         surf(xpi,ypi)%dem = newdem
      end do
   end do


   return
end subroutine complex_wall_texture


subroutine complex_terrace(user,surf,crater,deltaMtot)
   !Makes terraced walls by applying noisy topographic diffusion in discrete radial zones along the wall
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),intent(inout) :: deltaMtot

   ! Internal variables
   real(DP) :: newdem,elchange,rad,xbar,ybar,r,flr,hprof,tprof
   integer(I4B) :: xpi,ypi,i,j,inc,octave
   real(DP), dimension(2) :: rn

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   real(DP) :: noise,dnoise,tnoise

   ! Complex crater wall
   integer(I4B)            :: nscallops             ! Approximate number of scallop features on edge of crater
   integer(I4B),parameter  :: offset = 3000         ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP)                :: scallop_width          ! Vertical height of noise features as a function of crater radius at the first octave
                                                    ! Higher values make the scallop features broader
   real(DP)                :: scallop_p           ! The slope of the power law function that defines the noise shape for scallop features
                                                    ! Higher values makes the inner sides of the scallop features more rounded, lower values 
                                                    ! make them have sharper points at the inflections
   integer(I4B)                  :: nterraces       ! Number of terraces
   real(DP)                      :: terracefac      ! Power law scaling for relative terrace sizes
   integer(I4B)                  :: terrace         ! The current terrace
   real(DP)                      :: router          ! The radius of the terrace outer edge
   real(DP)                      :: rinner          ! The radius of the terrace outer edge
   real(DP)                      :: upshift,dfloor
   real(DP)                      :: h_scallop_profile ! Profile of the scallop wall
   integer(I4B)                  :: num_oct_tfloor 
   real(DP)                      :: noise_height_tfloor
   real(DP)                      :: freq_tfloor
   real(DP)                      :: pers_tfloor
   real(DP)                      :: xy_size_tfloor
   real(DP)                      :: rimfloor
   logical                       :: isterrace
   
   !Executable code
   call random_number(rn)


   
   ! Copernicus values
   terracefac = 3.0_DP
   nterraces = 8
   nscallops = 16 
   scallop_p = 0.6_DP 
   scallop_width = 0.20_DP 
   h_scallop_profile = 2.0_DP
   rimfloor = crater%melev + 0.5_DP * crater%rimheight
   upshift = 0.0_DP

   num_oct_tfloor  = 4
   noise_height_tfloor = crater%floordepth / nterraces
   freq_tfloor = 1.5_DP
   pers_tfloor = 0.5_DP
   xy_size_tfloor = 5.0_DP / crater%fcrat 

   ! Lansberg values
   !terracefac = 1.0_DP
   !nterraces = 16
   !nscallops = 8  
   !scallop_p = 0.6_DP
   !scallop_width = 0.10_DP 
   !^^^^^^^^^^^^^^^

   rad = 2.0_DP * crater%frad
 
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

            xbar = xpi * user%pix - crater%xl 
            ybar = ypi * user%pix - crater%yl

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            newdem = surf(xpi,ypi)%dem

            r = sqrt(xbar**2 + ybar**2) / crater%frad

            ! Make scalloped terraces
            znoise = scallop_width**(1._DP / (2 * scallop_p)) 
            xynoise = (nscallops / PI) / crater%fcrat
            dnoise = util_perlin_noise(xynoise * xbar + terrace * offset * rn(1), &
                                       xynoise * ybar + terrace * offset * rn(2)) * znoise
            noise = (dnoise**2)**scallop_p
            hprof = (r / router)**(-1)

            ! Make textured floor of terrace 
            tnoise = 0.0_DP
            do octave = 1, num_oct_tfloor
               xynoise = xy_size_tfloor * freq_tfloor ** (octave - 1) 
               znoise = noise_height_tfloor  * (pers_tfloor ) ** (octave - 1) 
               tnoise = tnoise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                                 xynoise * ybar + offset * rn(2))* znoise
            end do

            isterrace = 1.0_DP - noise < hprof

            if (r < 1.0_DP) then
               tprof = crater_profile(user,crater,rinner) * (r/rinner)**h_scallop_profile + tnoise ! This is the floor profile that replaces the old one at each terrace
            else
               tprof = (crater_profile(user,crater,rinner) + crater%ejrim * (rinner)**(-EJPROFILE)) + tnoise ! This is the floor profile that replaces the old one at each terrace
            end if

            if (isterrace) then
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

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         r = sqrt(xbar**2 + ybar**2) / crater%frad

         ! Make scalloped rim
         znoise = (scallop_width)**(1._DP / (2 * scallop_p))
         xynoise = (nscallops / PI) / crater%fcrat
         dnoise = util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                    xynoise * ybar + offset * rn(2)) * znoise
         noise = (dnoise**2)**scallop_p

         hprof = r**(-1)
         isterrace = 1.0_DP - noise < hprof

         ! Make textured floor of rim
         tnoise = 0.0_DP
         do octave = 1, num_oct_tfloor
            xynoise = xy_size_tfloor * freq_tfloor ** (octave - 1) 
            znoise = noise_height_tfloor  * (pers_tfloor ) ** (octave - 1) 
            tnoise = tnoise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2))* znoise
         end do

         tprof = rimfloor + tnoise

         if (isterrace) then
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

         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         newdem = surf(xpi,ypi)%dem

         r = sqrt(xbar**2 + ybar**2) / crater%frad
         
         if (r > flr) then
            newdem = newdem + upshift * max(min(3.0_DP - 2 * r**2,1.0_DP),0.0_DP)
            elchange  = newdem - surf(xpi,ypi)%dem
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         end if
      end do 
   end do

   return
end subroutine complex_terrace


subroutine ejecta_texture(user,surf,crater,deltaMtot,inc,ejecta_dem)
   ! Adds realistic texture to the ejecta
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),intent(in) :: deltaMtot
   integer(I4B),intent(in) :: inc
   real(DP),dimension(-inc:inc,-inc:inc),intent(inout) :: ejecta_dem

   ! Internal variables
   real(DP) :: newdem,elchange,xbar,ybar,xstretch,ystretch,ejbmass,ejbmassnew,fmasscons,r,phi,areafrac
   integer(I4B) :: i,j,xpi,ypi
   integer(I4B),dimension(2,-inc:inc,-inc:inc) :: indarray
   real(DP), dimension(2) :: rn

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise, hprof,xysplat,zsplat,dsplat,splatnoise
   integer(I4B) :: octave
   real(DP) :: noise,dnoise

   ! Ejecta base texture parameters
   integer(I4B)        :: num_octaves        ! Number of Perlin noise octaves
   integer(I4B)        :: offset        ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP)            :: xy_noise_fac           ! Spatial "size" of noise features at the first octave
   real(DP)            :: noise_height             ! Spatial "size" of noise features at the first octave
   real(DP)            :: freq              ! Spatial size scale factor multiplier at each octave level
   real(DP)            :: pers            ! The relative size scaling at each octave level


   ! Splat pattern
   integer(I4B)            :: nsplats             ! Approximate number of scallop features on edge of crater
   real(DP)                :: splat_height            ! Vertical height of noise features as a function of crater radius at the first octave
                                                    ! Higher values make the splat features broader
   real(DP)                :: splat_slope           ! The slope of the power law function that defines the noise shape for splat features
                                                    ! Higher values makes the inner sides of the splat features more rounded, lower values 
                                                    ! make them have sharper points at the inflections
   real(DP)                :: splat_stretch           ! Power law factor that stretches the splat pattern out
                                                     ! Higher values make the splat more stretched
   logical                       :: insplat
   real(DP)                      :: splatmag           ! The magnitude of the splat features relative to the ejecta thickness
   integer(I4B)         :: nsplat_octaves

   !Executable code
   call random_number(rn)


   ! Copernicus values
   num_octaves = 5
   offset = 7800
   xy_noise_fac = 6.0_DP
   noise_height = 4.0_DP
   freq = 2.0_DP
   pers = 0.5_DP

   nsplats = 32
   nsplat_octaves = 4
   splat_height = 14.0_DP
   splat_slope = 0.8_DP
   splat_stretch = 16.0_DP
   splatmag = 0.10_DP
   
   ! open(unit=12,file='params.txt',status='old')
   ! read(12,*) num_octaves
   ! read(12,*) xy_noise_fac
   ! read(12,*) noise_height
   ! close(12)

   ! Get the ejecta mass
   ejbmass = sum(ejecta_dem)
   if (ejbmass <= VSMALL) return

   ! First strip away the original ejecta from the surface

   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         surf(xpi,ypi)%dem = surf(xpi,ypi)%dem - ejecta_dem(i,j)

         ! Save index map
         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi
      end do
   end do


   ! Add the base texture to the ejecta proportional to the thickness
   do j = -inc,inc
      do i = -inc,inc

         xpi = indarray(1,i,j)
         ypi = indarray(2,i,j)

         xbar = (crater%xlpx + i) * user%pix - crater%xl 
         ybar = (crater%ylpx + j) * user%pix - crater%yl

         r = sqrt(xbar**2 + ybar**2) / crater%frad
         phi = atan2(ybar,xbar)

         ! Apply stretch to splats
         xstretch = r**(1.0_DP / splat_stretch) * cos(phi) * crater%frad
         ystretch = r**(1.0_DP / splat_stretch) * sin(phi) * crater%frad
  
         areafrac = util_area_intersection(inc * user%pix,xbar,ybar,user%pix) 
         areafrac = areafrac * (1.0_DP - util_area_intersection(crater%frad,xbar,ybar,user%pix))
      
         areafrac = areafrac * (1.0_DP - max(min(2._DP - r,1.0_DP),0.0_DP)) ! Blend in with the wall texture

         ! Make the splat pattern
         splatnoise = 0.0_DP
         do octave = 1,nsplat_octaves
            xysplat = (nsplats / PI) * freq ** (octave -1) / crater%fcrat 
            zsplat= splat_height**(1._DP / (2 * splat_slope)) * (pers ) ** (octave - 1)
            dsplat = util_perlin_noise(xysplat * xstretch + offset * rn(1), &
                                       xysplat * ystretch + offset * rn(2)) * zsplat
            dsplat = (dsplat**2)**(splat_slope)
            splatnoise = splatnoise + dsplat
         end do
         hprof =  r**(1.0_DP)
         insplat = 1.0_DP + splatnoise > hprof

         ! make base texture and then add extra layers if we are in one
         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) / crater%fcrat 
            znoise = noise_height  * (pers ) ** (octave - 1) * ejecta_dem(i,j)  
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2))* znoise

            if (insplat) noise = noise + (1.0_DP + splatnoise - hprof) * ejecta_dem(i,j) * splatmag
         end do

         ejecta_dem(i,j) = max(ejecta_dem(i,j) + noise * areafrac,0.0_DP)

      end do
   end do

   ! Save new total mass for mass conservation calculation
   ejbmassnew = sum(ejecta_dem)

   ! Rescale ejecta blanket to conserve mass
   fmasscons = (ejbmass - deltaMtot) / ejbmassnew
   ejecta_dem(:,:) = ejecta_dem(:,:) * fmasscons

   !Put the ejecta back onto the surface
   do j = -inc,inc
      do i = -inc,inc

         xpi = indarray(1,i,j)
         ypi = indarray(2,i,j)

         surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + ejecta_dem(i,j)

      end do
   end do

   return
end subroutine ejecta_texture


subroutine crater_realistic_slope_texture(user,critical_value,inc,critarray)
   ! Adds noise to the critical slope to give texture to regions that undergo slope collapse
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   real(DP),intent(in) :: critical_value
   integer(I4B),intent(in) :: inc
   real(DP),dimension(-inc:inc,-inc:inc),intent(out) :: critarray

   ! Internal variables
   real(DP) :: xynoise, znoise,xbar,ybar,r,noise
   integer(I4B) :: octave,i,j
   real(DP), dimension(2) :: rn

   ! Topographic noise parameters
   integer(I4B), parameter :: num_octaves  = 4   ! Number of Perlin noise octaves
   integer(I4B), parameter :: offset = 4000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: xy_noise_fac = 0.125_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: noise_height = 0.4e0_DP  ! Magnitude of noise features at the first octave
   real(DP), parameter :: freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: pers = 1.00_DP  ! The relative size scaling at each octave level

   ! Executable code
   call random_number(rn)

   do j = -inc,inc
      do i = -inc,inc

         xbar = real(i,kind=DP)
         ybar = real(j,kind=DP)

         r = sqrt(xbar**2 + ybar**2) 

         noise = 0.0_DP
         do octave = 1, num_octaves 
            xynoise = xy_noise_fac * freq ** (octave - 1) 
            znoise  = noise_height * pers ** (octave - 1) 
            noise = noise + util_perlin_noise(xynoise * xbar + offset * rn(1), &
                                              xynoise * ybar + offset * rn(2)) * znoise
         end do
         
         critarray(i,j) = max(critical_value * (1.0_DP + noise),0.0_DP)
      end do
   end do

   return
end subroutine crater_realistic_slope_texture

