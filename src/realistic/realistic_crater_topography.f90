!***** realistic/realistic_crater_topography
! Name
!   realistic_crater_topography 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_crater
!   * module_realistic
!   
!   call realistic_crater_topography(user,surf,crater,domain,deltaMtot)
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
subroutine realistic_crater_topography(user,surf,crater,domain,ejecta_dem)
   use module_globals
   use module_util
   use module_crater
   use module_realistic!, EXCEPT_THIS_ONE => realistic_crater_topography
   implicit none

   ! in and out
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(inout) :: ejecta_dem

   ! internal 
   real(DP) :: deltaMtot
   integer(I4B) :: inc
 
   ! excutable
   deltaMtot = 0.0_DP
   call realistic_rim(user,surf,crater,deltaMtot)
   
   return
end subroutine realistic_crater_topography

 

 

! new subroutines added by jundu on 9/28/2022




subroutine realistic_rim(user,surf,crater,deltaMtot)
    
   use module_globals
   use module_util
   use module_crater
   use module_realistic!, EXCEPT_THIS_ONE => realistic_rim
   implicit none

   ! in and out
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),intent(inout) :: deltaMtot
   
   ! internal
   type(psdtype) :: psd_distance,psd_elevation                                                  ! PSDs of rim distance and elevation
   real(DP) :: newdem,elchange,rad_roi,xbar,ybar
   integer(I4B) :: xpi,ypi,i,j,inc 
   real(DP) :: arc_length,radial_distance,random
   real(DP) :: rim_distance,rim_distance_delta,rim_elevation,rim_elevation_delta                ! rim distance and elevation and their variations. PSD only characterizes the variation.
   real(DP) :: diameter_in_m,gridsize
   real(DP) :: rim_elevation_old 
   integer(I4B),parameter :: infile = 14
   integer(I4B),parameter :: outfile = 15
   real(DP),dimension(:), allocatable :: amplitude_distance,wavelength_distance,phase_distance
   real(DP),dimension(:), allocatable :: amplitude_elevation,wavelength_elevation,phase_elevation
   integer(I4B),dimension(:), allocatable ::random_seed_array
   integer(I4B) :: n
    
   ! executable  
   gridsize=user%pix
   diameter_in_m=crater%fcrat
   rad_roi = 2.0_DP * crater%frad 
   inc = max(min(nint(rad_roi / user%pix),PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)

   ! PSDs of rim distance
   psd_distance%diameter_in_km=diameter_in_m/1000                                        
   psd_distance%num_vertices=5000 
   psd_distance%num_sine=psd_distance%num_vertices/2
   psd_distance%input_x_max=PI*psd_distance%diameter_in_km
   psd_distance%diameter_in_km_trans=20.0_DP
   psd_distance%bp2_x_k_b_sigma= [0.06294416243654823_DP , -0.23388324873096447_DP, 0.0056111111111111145_DP, 0.9127777777777776_DP, 0.06082223144489516_DP]
   psd_distance%bp2_y_k_b_sigma= [0.17451428571428573_DP , 0.5097142857142856_DP  , 0.0015652173913043492_DP, 3.968695652173913_DP , 0.39645422747641196_DP]
   psd_distance%bp3_y_k_b_sigma= [0.13316666666666666_DP , 1.836666666666667_DP   , 0.0004552129221732701_DP, 4.4908957415565345_DP, 0.27156105098740746_DP]
   psd_distance%bp4_y_k_b_sigma= [0.200561797752809_DP   , -2.4412359550561797_DP , 0.005084745762711865_DP , 1.4683050847457628_DP, 0.349619629586608_DP  ]
   psd_distance%slope1_k_b_sigma=[0.004916013437849945_DP, 2.9333146696528556_DP  , -99999.0_DP             , -99999.0_DP          , 0.36717269233927413_DP]
   psd_distance%misfit_k_b_sigma(1)=0.000000_DP

   ! PSDs of rim elevation
   psd_elevation%diameter_in_km=diameter_in_m/1000.0_DP                                       
   psd_elevation%num_vertices=5000 
   psd_elevation%num_sine=psd_elevation%num_vertices/2
   psd_elevation%input_x_max=PI*psd_elevation%diameter_in_km
   psd_elevation%diameter_in_km_trans=20.0_DP
   psd_elevation%bp2_x_k_b_sigma=[0.05533333333333333_DP, -0.10666666666666666_DP, 0.0040595399188092015_DP, 0.918809201623816_DP,0.11430090458471147_DP]
   psd_elevation%bp2_y_k_b_sigma=[0.12594736842105264_DP, -0.2859473684210526_DP, 0.003978349120433018_DP, 2.1534330175913396_DP, 0.40875902323519153_DP]
   psd_elevation%bp3_y_k_b_sigma=[0.0_DP, 4.0_DP, 0.0_DP, 4.0_DP,0.7435034294351548_DP]
   psd_elevation%bp4_y_k_b_sigma(1)=-99999.0_DP
   psd_elevation%slope1_k_b_sigma=[0.09919786096256683_DP, 1.7070427807486632_DP, -0.008135593220338983_DP, 3.8537118644067796_DP,0.42722290048928596_DP]
   psd_elevation%misfit_k_b_sigma(1)=0.000000_DP



 
   call Calculate_am_wl_phase_from_diameter(psd_distance,amplitude_distance,wavelength_distance,phase_distance)
   call Calculate_am_wl_phase_from_diameter(psd_elevation,amplitude_elevation,wavelength_elevation,phase_elevation)




   rim_elevation_old=crater%rimheight 
 
   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j
         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         arc_length=(atan2(real(j,8),real(i,8))+PI)*diameter_in_m/2.0_DP/1000.0_DP
         radial_distance=  sqrt(xbar**2 + ybar**2) 

         call Create_rim(arc_length,psd_distance,amplitude_distance,wavelength_distance,phase_distance,rim_distance_delta)
         call Create_rim(arc_length,psd_elevation,amplitude_elevation,wavelength_elevation,phase_elevation,rim_elevation_delta)
         
         rim_distance=rim_distance_delta+diameter_in_m/2.0_DP
         rim_elevation=rim_elevation_delta+rim_elevation_old
         crater%rimheight=rim_elevation     

         if (radial_distance>rim_distance) then 
            newdem=crater_profile(user,crater,radial_distance/rim_distance)
            elchange  = newdem-surf(xpi,ypi)%dem       
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         endif

         if (radial_distance<rim_distance .and. radial_distance>crater%floordiam/2.0_DP) then 
            newdem=crater_profile(user,crater,radial_distance/rim_distance)
            elchange  = newdem-surf(xpi,ypi)%dem       
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         endif

      end do
   end do


    

end subroutine realistic_rim



 ! Universal shapde model from PSD

subroutine Calculate_am_wl_phase_from_diameter(psd_1D,amplitude,wavelength,phase)
   use module_globals
   use module_realistic
   implicit none
   ! in and out
   type(psdtype),intent(inout) :: psd_1D
   real(DP),dimension(:), allocatable,intent(out) :: amplitude,wavelength,phase
   ! internal  
   integer(I4B) :: i
   real(DP),dimension(:), allocatable :: psd
   integer(I4B),parameter :: infile = 14
   integer(I4B),parameter :: outfile = 15

   ! excutable
   call Calculate_breakpoint_slope_from_diameter(psd_1D)
   call Calculate_targetPSD_from_breakpoint_slope(psd_1D,wavelength,psd)
   call Calculate_am_wl_phase_from_targetPSD(psd_1D,wavelength,psd,amplitude,phase)

   open(outfile,file="output_psd.txt",status='replace')
      do i = 1,  size(wavelength)
         write(outfile,'(ES18.10,1(",",ES18.10))') wavelength(i),psd(i)
      end do
   close(outfile)

end subroutine Calculate_am_wl_phase_from_diameter

subroutine Calculate_breakpoint_slope_from_diameter(psd_1D)
   use module_globals
   use module_realistic
   use module_util
   implicit none
   ! in and out
   type(psdtype),intent(inout) :: psd_1D
   ! internal
   real(DP) :: random_number_normal
   
   ! excutable
   if (  psd_1D%diameter_in_km  <  psd_1D%diameter_in_km_trans  ) then
      psd_1D%bp2_x=psd_1D%bp2_x_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp2_x_k_b_sigma(2)  
      psd_1D%bp2_y=psd_1D%bp2_y_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp2_y_k_b_sigma(2)
      psd_1D%bp3_y=psd_1D%bp3_y_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp3_y_k_b_sigma(2)
      psd_1D%bp4_y=psd_1D%bp4_y_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp4_y_k_b_sigma(2)
   
   else

      psd_1D%bp2_x=psd_1D%bp2_x_k_b_sigma(3)*psd_1D%diameter_in_km+psd_1D%bp2_x_k_b_sigma(4)
      psd_1D%bp2_y=psd_1D%bp2_y_k_b_sigma(3)*psd_1D%diameter_in_km+psd_1D%bp2_y_k_b_sigma(4)
      psd_1D%bp3_y=psd_1D%bp3_y_k_b_sigma(3)*psd_1D%diameter_in_km+psd_1D%bp3_y_k_b_sigma(4)
      psd_1D%bp4_y=psd_1D%bp4_y_k_b_sigma(3)*psd_1D%diameter_in_km+psd_1D%bp4_y_k_b_sigma(4)

   endif


   if   (psd_1D%slope1_k_b_sigma(3)>-100.0_DP) then
      if (  psd_1D%diameter_in_km  <  psd_1D%diameter_in_km_trans  ) then
         psd_1D%slope1=psd_1D%slope1_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%slope1_k_b_sigma(2)
      else
         psd_1D%slope1=psd_1D%slope1_k_b_sigma(3)*psd_1D%diameter_in_km+psd_1D%slope1_k_b_sigma(4)
      endif
   else
      psd_1D%slope1=psd_1D%slope1_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%slope1_k_b_sigma(2)
   endif


   psd_1D%bp3_y=7
   print *,'bp3_y: ',psd_1D%bp3_y

   ! call util_random_number_normal(random_number_normal)
   ! print *,'1',random_number_normal
   ! psd_1D%bp2_x= psd_1D%bp2_x  + random_number_normal * psd_1D%bp2_x_k_b_sigma(5)   
   ! call util_random_number_normal(random_number_normal)
   ! print *,'2',random_number_normal
   ! psd_1D%bp2_y= psd_1D%bp2_y  + random_number_normal * psd_1D%bp2_y_k_b_sigma(5) 
   ! call util_random_number_normal(random_number_normal)
   ! print *,'3',random_number_normal
   ! psd_1D%bp3_y= psd_1D%bp3_y  + random_number_normal * psd_1D%bp3_y_k_b_sigma(5) 
   ! call util_random_number_normal(random_number_normal)
   ! print *,'4',random_number_normal
   ! psd_1D%bp4_y= psd_1D%bp4_y  + random_number_normal * psd_1D%bp4_y_k_b_sigma(5) 
   ! call util_random_number_normal(random_number_normal)
   ! print *,'5',random_number_normal
   ! psd_1D%slope1=psd_1D%slope1 + random_number_normal * psd_1D%slope1_k_b_sigma(5) 

end subroutine Calculate_breakpoint_slope_from_diameter

subroutine Calculate_targetPSD_from_breakpoint_slope(psd_1D,wavelength,psd)
   use module_globals
   use module_realistic
   use module_util
   implicit none
   !in and out
   type(psdtype), intent(in)  :: psd_1D
   real(DP) ,dimension(:),allocatable,intent(out) :: wavelength,psd 
   ! internal
   integer(I4B) ::  i,bp2_x_index,bp3_x_index
   real(DP) ::  bp3_x,bp4_x
   real(DP) ::  slope_12,intercept_12,slope_23,intercept_23,slope_34,intercept_34
   real(DP) :: random_number_normal

   
   ! excutable
   allocate(wavelength(psd_1D%num_sine))
   allocate(psd(psd_1D%num_sine))
  !----------------------------------------------------------------------------------------------------------------------
   do i = 1, psd_1D%num_sine
      wavelength(i) =  psd_1D%input_x_max / (psd_1D%num_sine-i+1)
   end do
  !----------------------------------------------------------------------------------------------------------------------
      
   if (psd_1D%bp4_y_k_b_sigma(1)>-100.0_DP  ) then 
   !-----------------------------------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------------------------------------------
      ! 3 segments, 4 breakpoints
      bp3_x=log10(psd_1D%input_x_max/2.0_DP) 
      bp4_x=log10(psd_1D%input_x_max) 

      do i = 1, psd_1D%num_sine
         bp2_x_index = i
         if (wavelength(i) > 10**psd_1D%bp2_x) exit
      end do
      do i = 1, psd_1D%num_sine
         bp3_x_index = i
         if (wavelength(i) > 10**bp3_x) exit
      end do
   !----------------------------------------------------------------------------------------------------------------------
      slope_12=psd_1D%slope1
      intercept_12 = psd_1D%bp2_y - slope_12 * psd_1D%bp2_x 

      slope_23=(psd_1D%bp3_y-psd_1D%bp2_y)/(bp3_x-psd_1D%bp2_x)
      intercept_23 = psd_1D%bp3_y - slope_23 * bp3_x 

      slope_34=(psd_1D%bp4_y-psd_1D%bp3_y)/(bp4_x-bp3_x)
      intercept_34 = psd_1D%bp4_y - slope_34 * bp4_x 
   !----------------------------------------------------------------------------------------------------------------------
      do i = 1, bp2_x_index-1
         call util_random_number_normal(random_number_normal)
         psd(i) = 10** (   slope_12*log10(wavelength(i))    +  intercept_12  )!+  random_number_normal*psd_1D%misfit_k_b_sigma(1)  )      
      end do     
      do i = bp2_x_index, bp3_x_index-1
         call util_random_number_normal(random_number_normal)
         psd(i) = 10** (   slope_23*log10(wavelength(i))    +  intercept_23 )! +  random_number_normal*psd_1D%misfit_k_b_sigma(1)  )
      end do     
      do i = bp3_x_index, psd_1D%num_sine
         call util_random_number_normal(random_number_normal)
         psd(i) = 10** (   slope_34*log10(wavelength(i))    +  intercept_34 )! +  random_number_normal*psd_1D%misfit_k_b_sigma(1)  )
      end do

   else
   !-----------------------------------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------------------------------------------
      ! 2 segments, 3 breakpoints
      bp3_x=log10(psd_1D%input_x_max) 
 
      do i = 1, psd_1D%num_sine
         bp2_x_index = i
         if (wavelength(i) > 10**psd_1D%bp2_x) exit
      end do
   !----------------------------------------------------------------------------------------------------------------------
      slope_12=psd_1D%slope1
      intercept_12 = psd_1D%bp2_y - slope_12 * psd_1D%bp2_x 

      slope_23=(psd_1D%bp3_y-psd_1D%bp2_y)/(bp3_x-psd_1D%bp2_x)
      intercept_23 = psd_1D%bp3_y - slope_23 * bp3_x 
   !----------------------------------------------------------------------------------------------------------------------
      do i = 1, bp2_x_index-1
         call util_random_number_normal(random_number_normal)
         psd(i) = 10** (   slope_12*log10(wavelength(i))    +  intercept_12  )!+  random_number_normal*psd_1D%misfit_k_b_sigma(1)  )
      end do     
      do i = bp2_x_index, psd_1D%num_sine
         call util_random_number_normal(random_number_normal)
         psd(i) = 10** (   slope_23*log10(wavelength(i))    +  intercept_23 )! +  random_number_normal*psd_1D%misfit_k_b_sigma(1)  )
      end do
   endif
end subroutine Calculate_targetPSD_from_breakpoint_slope

subroutine Calculate_am_wl_phase_from_targetPSD(psd_1D,wavelength,psd,amplitude,phase)
   use module_globals
   use module_realistic
   implicit none
   !in and out
   type(psdtype), intent(in)  :: psd_1D
   real(DP) ,dimension(:),intent(in) :: wavelength,psd 
   real(DP) ,dimension(:),allocatable,intent(out) :: amplitude,phase 
   ! internal
   integer(I4B) ::  i
   real(DP) ::  phase_random
   ! excutable
   allocate(amplitude(psd_1D%num_sine))
   allocate(phase(psd_1D%num_sine))

  !----------------------------------------------------------------------------------------------------------------------
   do i = 1, psd_1D%num_sine
      amplitude(i)=  sqrt(psd(i) * psd_1D%input_x_max / (psd_1D%num_vertices ** 2))
      call RANDOM_NUMBER(phase_random)
      phase(i)=wavelength(i) *phase_random   
   end do
   
end subroutine Calculate_am_wl_phase_from_targetPSD

subroutine Create_rim(arc_length,psd_1D,amplitude,wavelength,phase,rim_parameter)
   use module_globals
   use module_realistic
   implicit none
   ! in and out
   real(DP),intent(in) :: arc_length
   type(psdtype),intent(in) :: psd_1D
   real(DP),dimension(:),intent(in) :: amplitude
   real(DP),dimension(:),intent(in) :: wavelength
   real(DP),dimension(:),intent(in) :: phase
   real(DP),intent(out) :: rim_parameter
   ! internal
   integer(I4B) :: i
   real(DP) :: rim_parameter_ind

   ! excutable
   rim_parameter = 0.0_DP
   do i = 1, psd_1D%num_sine
      rim_parameter_ind= amplitude(i) * sin( 2 * PI* 1/wavelength(i) * (  arc_length-   phase(i)   )     )      
      rim_parameter=rim_parameter+rim_parameter_ind
   end do
   rim_parameter=rim_parameter*1000.0_DP

end subroutine Create_rim
  