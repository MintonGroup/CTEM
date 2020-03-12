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
   real(DP),intent(out) :: deltaMtot

   ! Internal variables
   real(DP) :: lradsq,newelev, x_relative, y_relative 
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq,deltaMi,rimheight
   logical :: lastloop
   real(DP), dimension(2) :: rn


   ! Topographic noise parameters
   real(DP) :: xynoise, znoise
   integer(I4B) :: octave
   real(DP) :: noise


   ! Complex crater floor
   integer(I4B), parameter :: complex_floor_num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: complex_floor_offset = 10000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: complex_floor_xy_noise_fac = 2.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: complex_floor_noise_height = 1e-3_DP ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP), parameter :: complex_floor_freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: complex_floor_pers = 0.80_DP  ! The relative size scaling at each octave level

   ! Complex crater peak
   integer(I4B), parameter :: complex_peak_num_octaves  = 8   ! Number of Perlin noise octaves
   integer(I4B), parameter :: complex_peak_offset = 1000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: complex_peak_rad = 0.20_DP      ! Radial size of central peak relative to frad
   real(DP), parameter :: complex_peak_xy_noise_fac = 16.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: complex_peak_freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: complex_peak_pers = 0.80_DP  ! The relative size scaling at each octave level

   ! Complex crater wall
   integer(I4B), parameter :: complex_wall_num_octaves  = 3   ! Number of Perlin noise octaves
   integer(I4B), parameter :: complex_wall_offset = 2000 ! Scales the random xy-offset so that each crater's random noise is unique 
   real(DP), parameter :: complex_wall_xy_noise_fac = 2.0_DP  ! Spatial "size" of noise features at the first octave
   real(DP), parameter :: complex_wall_freq = 2.0_DP     ! Spatial size scale factor multiplier at each octave level
   real(DP), parameter :: complex_wall_pers = 0.80_DP  ! The relative size scaling at each octave level
   real(DP), parameter :: complex_wall_noise_height = 0.10_DP ! Vertical height of noise features as a function of crater radius at the first octave
   real(DP), parameter :: complex_wall_slope = 25._DP * DEG2RAD
   integer(I4B), parameter :: complex_wall_terracefac = 16 ! Spatial size factor for terraces relative to crater radius
   integer(I4B) :: complex_wall_terrace_num


   ! Executable code

   call random_number(rn)

   ! determine area to effect
   ! First make the interior of the crater
   inc = max(min(crater%fradpx,PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   deltaMtot = 0.0_DP 
   incsq = inc**2
   
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         iradsq = i**2 + j**2 
         ! find distance from crater center
         
         ! find elevation and grid point
         newelev = crater%melev + ((i * crater%xslp) + (j * crater%yslp)) * user%pix
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




         deltaMtot = deltaMtot + deltaMi

      end do
   end do


!            ! Make the floor topographically "noisy"
!            noise = 0.0_DP
!            do octave = 1, complex_floor_num_octaves
!               xynoise = complex_floor_xy_noise_fac * complex_floor_freq ** octave / crater%fcrat
!               znoise = complex_floor_noise_height * complex_floor_pers ** (octave - 1) * crater%fcrat
!               noise = noise + util_perlin_noise(xynoise * x_relative + complex_floor_offset * rn(1), &
!                                                 xynoise * y_relative + complex_floor_offset * rn(2)) * znoise
!            end do
!            newdem = newdem + max(noise,0.0_DP)
!
!            ! Add in "noisy" central peak
!            if (r < complex_peak_rad) then
!               Hpeak = 0.032e3_DP * (crater%fcrat * 1e-3_DP)**(0.900_DP) ! Pike (1977)
!               noise = 0.0_DP
!               do octave = 1, complex_peak_num_octaves 
!                  xynoise = complex_peak_xy_noise_fac * complex_peak_freq ** octave / crater%fcrat
!                  znoise = Hpeak * (1.0_DP - r / complex_peak_rad) * complex_peak_pers ** (octave - 1)
!                  noise = noise + util_perlin_noise(xynoise * x_relative + complex_peak_offset * rn(1), &
!                                                    xynoise * y_relative + complex_peak_offset * rn(2)) * znoise
!               end do 
!               newdem = newdem + max(noise,0.0_DP) 
!            end if
!         else ! Terraced walls 
!            complex_wall_terrace_num = int(complex_wall_terracefac * r) + 1
!            noise = 0.0_DP
!            do octave = 1, complex_wall_num_octaves
!               xynoise = complex_wall_xy_noise_fac * complex_wall_freq ** octave / crater%fcrat
!               znoise = complex_wall_noise_height * complex_wall_pers ** (octave - 1) * crater%fcrat
!               noise = noise + util_perlin_noise(xynoise * x_relative + complex_wall_terrace_num * complex_wall_offset * rn(1), &
!                                                 xynoise * y_relative + complex_wall_terrace_num * complex_wall_offset * rn(2)) &
!                                                 * znoise
!            end do
!            newdem = max(newdem + noise,crater%melev - pikeD) !/ cos(complex_wall_slope)


   return
end subroutine crater_realistic_topography

