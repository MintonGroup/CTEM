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
   real(DP) :: lradsq,x_relative, y_relative 
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   real(DP) :: r,xp,yp,fradsq,deltaMi,rimheight
   logical :: lastloop
   real(DP), dimension(2) :: rn

   interface
      subroutine complex_floor(user,surfi,crater,r,x_relative,y_relative,rn,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: r,x_relative, y_relative 
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(out) :: deltaMi
      end subroutine complex_floor

      subroutine complex_peak(user,surfi,crater,r,x_relative,y_relative,rn,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: r,x_relative, y_relative 
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(out) :: deltaMi
      end subroutine complex_peak

      subroutine complex_wall(user,surfi,crater,r,x_relative,y_relative,rn,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: r,x_relative, y_relative 
      real(DP),dimension(:),intent(in) :: rn
      real(DP),intent(out) :: deltaMi
      end subroutine complex_wall
   end interface

   ! Executable code
   call random_number(rn)

   ! determine area to effect
   ! First make the interior of the crater
   inc = max(min(crater%fradpx,PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   incsq = inc**2
   
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         iradsq = i**2 + j**2 
         ! find distance from crater center
         
         ! find elevation and grid point
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         xp = xpi * user%pix
         yp = ypi * user%pix
         
         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         x_relative = (crater%xl - xp)
         y_relative = (crater%yl - yp)
         
         lradsq = x_relative**2 + y_relative**2

         if (lradsq > crater%frad**2) cycle
         
         r = sqrt(x_relative**2+y_relative**2) / crater%frad

         if (crater%fcrat > crater%cxtran * 2) then
            call complex_floor(user,surf(xpi,ypi),crater,r,x_relative,y_relative,rn,deltaMi)
            call complex_peak (user,surf(xpi,ypi),crater,r,x_relative,y_relative,rn,deltaMi)
            call complex_wall (user,surf(xpi,ypi),crater,r,x_relative,y_relative,rn,deltaMi)
         end if

         deltaMtot = deltaMtot + deltaMi

      end do
   end do


   return
end subroutine crater_realistic_topography

subroutine complex_peak(user,surfi,crater,r,x_relative,y_relative,rn,deltaMi)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: r,x_relative, y_relative 
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(out) :: deltaMi

   ! Complex crater peak
   real(DP), parameter :: complex_peak_rad = 0.20_DP      ! Radial size of central peak relative to frad
   integer(I4B), parameter :: complex_peak_num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: complex_peak_offset = 1000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: complex_peak_xy_noise_fac = 16.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: complex_peak_freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: complex_peak_pers = 0.80_DP  ! The relative size scaling at each octave level

   ! Internal variables
   real(DP) :: newdem,elchange,pikeD,Hpeak

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise

   if (r > complex_peak_rad) return
   ! Add in "noisy" central peak
   newdem = surfi%dem
   Hpeak = 0.032e3_DP * (crater%fcrat * 1e-3_DP)**(0.900_DP) ! Pike (1977)
   noise = 0.0_DP
   do octave = 1, complex_peak_num_octaves 
      xynoise = complex_peak_xy_noise_fac * complex_peak_freq ** octave / crater%fcrat
      znoise = Hpeak * (1.0_DP - r / complex_peak_rad) * complex_peak_pers ** (octave - 1)
      noise = noise + util_perlin_noise(xynoise * x_relative + complex_peak_offset * rn(1), &
                                        xynoise * y_relative + complex_peak_offset * rn(2)) * znoise
   end do 
   newdem = newdem + max(noise,0.0_DP) 

   elchange  = newdem - surfi%dem
   deltaMi = elchange
   surfi%dem = newdem

end subroutine complex_peak

subroutine complex_floor(user,surfi,crater,r,x_relative,y_relative,rn,deltaMi)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: r,x_relative, y_relative 
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(out) :: deltaMi

   ! Complex crater floor parameters
   integer(I4B), parameter :: complex_floor_num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: complex_floor_offset = 10000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: complex_floor_xy_noise_fac = 2.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: complex_floor_noise_height = 1e-3_DP ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP), parameter :: complex_floor_freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: complex_floor_pers = 0.80_DP  ! The relative size scaling at each octave level

   ! Internal variables
   real(DP) :: newdem,elchange,pikeD,Hpeak

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise

   newdem = surfi%dem
   ! Make the floor topographically "noisy"
   noise = 0.0_DP
   do octave = 1, complex_floor_num_octaves
      xynoise = complex_floor_xy_noise_fac * complex_floor_freq ** octave / crater%fcrat
      znoise = complex_floor_noise_height * complex_floor_pers ** (octave - 1) * crater%fcrat
      noise = noise + util_perlin_noise(xynoise * x_relative + complex_floor_offset * rn(1), &
                                        xynoise * y_relative + complex_floor_offset * rn(2)) * znoise
   end do
   newdem = newdem + max(noise,0.0_DP)

   elchange  = newdem - surfi%dem
   deltaMi = elchange
   surfi%dem = newdem

end subroutine complex_floor


subroutine complex_wall(user,surfi,crater,r,x_relative,y_relative,rn,deltaMi)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_realistic_topography
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: r,x_relative, y_relative 
   real(DP),dimension(:),intent(in) :: rn
   real(DP),intent(out) :: deltaMi

   ! Internal variables
   real(DP) :: newdem,elchange,pikeD,pikeFloor,Hpeak

   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise

   ! Complex crater wall
   integer(I4B), parameter :: complex_wall_num_octaves  = 3   ! Number of Perlin noise octaves
   integer(I4B), parameter :: complex_wall_offset = 2000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: complex_wall_xy_noise_fac = 2.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: complex_wall_freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: complex_wall_pers = 0.80_DP  ! The relative size scaling at each octave level
   real(DP), parameter :: complex_wall_noise_height = 0.10_DP ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP), parameter :: complex_wall_slope = 25._DP * DEG2RAD
   integer(I4B), parameter :: complex_wall_terracefac = 32 ! Spatial size factor for terraces relative to crater radius
   integer(I4B) :: complex_wall_terrace_num


   newdem = surfi%dem
   pikeD = min(1.044e3_DP * (crater%fcrat * 1e-3_DP)**(0.301_DP), user%deplimit) ! Pike (1977) depth/diameter ratio of complex craters
   !pikeFloor = 
   pikeFloor = 0.7_DP
   if (r < pikeFloor) return
   complex_wall_terrace_num = int(complex_wall_terracefac * r) + 1
   noise = 0.0_DP
   do octave = 1, complex_wall_num_octaves
      xynoise = complex_wall_xy_noise_fac * complex_wall_freq ** octave / crater%fcrat
      znoise = complex_wall_noise_height * complex_wall_pers ** (octave - 1) * crater%fcrat
      noise = noise + util_perlin_noise(xynoise * x_relative + complex_wall_terrace_num * complex_wall_offset * rn(1), &
                                        xynoise * y_relative + complex_wall_terrace_num * complex_wall_offset * rn(2)) &
                                        * znoise
   end do
   newdem = max(newdem + noise,crater%melev - pikeD) !/ cos(complex_wall_slope)


   elchange  = newdem - surfi%dem
   deltaMi = elchange
   surfi%dem = newdem

end subroutine complex_wall




