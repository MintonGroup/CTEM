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
   use module_realistic, EXCEPT_THIS_ONE => realistic_crater_topography
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
   call make_realistic_crater(user,surf,crater,deltaMtot)
   
   return
end subroutine realistic_crater_topography

 

 

! new subroutines added by jundu on 9/28/2022




subroutine make_realistic_crater(user,surf,crater,deltaMtot)
    
   use module_globals
   use module_util
   use module_crater
   use module_realistic, EXCEPT_THIS_ONE => make_realistic_crater
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
   real(DP) :: floor_distance,floor_distance_delta
   real(DP) :: ejecta_distance,ejecta_distance_delta
   real(DP) :: diameter_in_m
   real(DP) :: rim_height_old,floor_diameter_old,rim_diameter_old
   integer(I4B),parameter :: infile = 14
   integer(I4B),parameter :: outfile = 15



   real(DP),dimension(:), allocatable :: amplitude_rim_distance,wavelength_rim_distance,phase_rim_distance
   real(DP),dimension(:), allocatable :: amplitude_floor_distance,wavelength_floor_distance,phase_floor_distance
   real(DP),dimension(:), allocatable :: amplitude_ejecta_distance,wavelength_ejecta_distance,phase_ejecta_distance
   real(DP),dimension(:), allocatable :: amplitude_rim_elevation,wavelength_rim_elevation,phase_rim_elevation
   integer(I4B),dimension(:), allocatable ::random_seed_array
   integer(I4B) :: n
    
   ! executable  

   diameter_in_m=crater%fcrat
   rad_roi = 2.0_DP * crater%frad 
   inc = max(min(nint(rad_roi / user%pix),PBCLIM*user%gridsize),1) + 1
    

   print *,user%pix,diameter_in_m,rad_roi,inc


   do j = 1,user%gridsize
      do i = 1,user%gridsize
         surf(j,i)%dem=0
      end do
   end do




   ! PSDs of rim distance
   psd_distance%diameter_in_km=diameter_in_m/1000                                        
   psd_distance%num_vertices=5000 
   psd_distance%num_sine=psd_distance%num_vertices/2
   psd_distance%input_x_max=PI*psd_distance%diameter_in_km
   psd_distance%diameter_in_km_trans=20.0_DP
   psd_distance%bp2_x_k_b_sigma= [0.036849225951533300_DP , 3.62699311558896E-25_DP, 0.38944510211326300_DP, 0.009322013219269230_DP, 0.5505442546452810_DP,0.3900022904707120_DP]
   psd_distance%bp2_y_k_b_sigma= [0.08090546246236460_DP , 0.9871469283188110_DP  , 0.5387009528903290_DP, 0.0199490078332091_DP , 2.20627602090192_DP,0.9702342937127700_DP]
   psd_distance%bp3_y_k_b_sigma= [0.0999302402137831_DP , 2.5221841961054200_DP   ,0.44782473707744800_DP, 0.003495544884342460_DP,4.450878102694230_DP,0.518635273427521_DP]
   psd_distance%bp4_y_k_b_sigma= [0.05303895069188420_DP   , 1.11858260986482E-20_DP , 1.092129604779800_DP , 0.017861817477899800_DP, 0.7035426642796880_DP,0.7265758599697800_DP  ]
   psd_distance%slope12_k_b_sigma=[-0.018167695145018800_DP, 3.1989973581963700_DP  , 0.5862827181492530_DP , -0.0006916915876435890_DP, 2.8494772870488700_DP,0.23144222028562200_DP]
   psd_distance%misfit_k_b_sigma=0.000000_DP
   call Calculate_am_wl_phase_from_diameter(psd_distance,amplitude_rim_distance,wavelength_rim_distance,phase_rim_distance)


   


   psd_distance%bp2_x_k_b_sigma= [0.0445647071_DP, 0.0249044164_DP, 0.2457300647_DP, 0.0033776036_DP, 0.8486464862_DP, 0.2592335367_DP]
   psd_distance%bp2_y_k_b_sigma= [0.1249485287_DP, 0.7508834243_DP, 0.8412738883_DP, 0.0018242070_DP, 3.2133698574_DP, 0.5984546470_DP]
   psd_distance%bp3_y_k_b_sigma= [0.1165535795_DP, 2.4968468812_DP, 0.3435521591_DP, 0.0007293887_DP, 4.8133306975_DP, 0.4967945371_DP]
   psd_distance%bp4_y_k_b_sigma= [0.1090009201_DP, 3.2535271225_DP, 0.8107845470_DP, 0.0032284355_DP, 5.3689768150_DP, 0.5578406988_DP]
   psd_distance%slope12_k_b_sigma=[-0.0333844825_DP, 3.5092973932_DP, 0.5485796877_DP, -0.0036202714_DP, 2.9140131695_DP, 0.2072551238_DP]
   psd_distance%misfit_k_b_sigma=0.000000_DP
   call Calculate_am_wl_phase_from_diameter(psd_distance,amplitude_floor_distance,wavelength_floor_distance,phase_floor_distance)



   psd_distance%bp2_x_k_b_sigma= [0.0209840178_DP, 0.0000000000_DP, 0.3137881706_DP, 0.0087635772_DP, 0.2444088123_DP, 0.2532110670_DP]
   psd_distance%bp2_y_k_b_sigma= [0.0196037572_DP, 3.5647808017_DP, 0.2978085502_DP, 0.0061716830_DP, 3.8334222875_DP, 0.5419519768_DP]
   psd_distance%bp3_y_k_b_sigma= [0.0592828357_DP, 4.7069303696_DP, 0.3348047038_DP, -0.0000033398_DP, 5.8926538800_DP, 0.2997192188_DP]
   psd_distance%bp4_y_k_b_sigma= [0.0615342686_DP, 5.0000000000_DP, 0.4785500634_DP, 0.0063998082_DP, 6.1026892077_DP, 0.6110155488_DP]
   psd_distance%slope12_k_b_sigma=[0.0073558947_DP, 2.6540309528_DP, 0.3107837181_DP, -0.0001997034_DP, 2.8051429155_DP, 0.4398559005_DP]
   psd_distance%misfit_k_b_sigma=0.000000_DP
   call Calculate_am_wl_phase_from_diameter(psd_distance,amplitude_ejecta_distance,wavelength_ejecta_distance,phase_ejecta_distance)





   ! PSDs of rim elevation
   psd_elevation%diameter_in_km=diameter_in_m/1000.0_DP                                       
   psd_elevation%num_vertices=5000 
   psd_elevation%num_sine=psd_elevation%num_vertices/2
   psd_elevation%input_x_max=PI*psd_elevation%diameter_in_km
   psd_elevation%diameter_in_km_trans=20.0_DP
   psd_elevation%bp2_x_k_b_sigma=[0.0417732075218209_DP , 2.48787783435609e-18_DP , 0.496231041304548_DP , 0.00210183162260986_DP , 0.793427517984220_DP , 0.308862248697524_DP]
   psd_elevation%bp2_y_k_b_sigma=[0.0593460236785858_DP , 0.408350328939867_DP , 1.53480069842064_DP , -0.00201958329120022_DP , 1.63566246833559_DP , 0.809757199553671_DP]
   psd_elevation%bp3_y_k_b_sigma=[0.0561324763037810_DP , 2.46101909781936_DP , 0.537503435665947_DP , -0.00292059238533839_DP , 3.64208047160175_DP , 0.462627466875698_DP]
   psd_elevation%bp4_y_k_b_sigma=[0.0159696581958859_DP , 3.25122353414224_DP , 0.660958440041299_DP , 0.00638368649794204_DP , 3.44294296810112_DP , 0.769783029462181_DP]
   psd_elevation%slope12_k_b_sigma=[0.0279866462789263_DP , 2.53280320608848_DP , 0.492438940196540_DP , 0.000799959919780047_DP , 3.07653693327141_DP , 0.400693610967074_DP]
   psd_elevation%misfit_k_b_sigma=0.000000_DP
   call Calculate_am_wl_phase_from_diameter(psd_elevation,amplitude_rim_elevation,wavelength_rim_elevation,phase_rim_elevation)




   rim_height_old=crater%rimheight 
   floor_diameter_old=crater%floordiam
   rim_diameter_old=crater%fcrat



 
   print *,crater%fcrat,crater%floordepth,crater%rimheight,crater%floordiam
   
   do j = -inc,inc
      do i = -inc,inc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j
         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         arc_length=(atan2(real(j,8),real(i,8))+PI)*diameter_in_m/2.0_DP/1000.0_DP
         radial_distance=  sqrt(xbar**2 + ybar**2) 

        
         call Make_outline(arc_length,psd_distance ,amplitude_rim_distance   ,wavelength_rim_distance   ,phase_rim_distance   ,rim_distance_delta)
         call Make_outline(arc_length,psd_distance ,amplitude_floor_distance ,wavelength_floor_distance ,phase_floor_distance ,floor_distance_delta)
         call Make_outline(arc_length,psd_distance ,amplitude_ejecta_distance,wavelength_ejecta_distance,phase_ejecta_distance,ejecta_distance_delta)
         call Make_outline(arc_length,psd_elevation,amplitude_rim_elevation  ,wavelength_rim_elevation  ,phase_rim_elevation  ,rim_elevation_delta)

         rim_distance=rim_diameter_old/2.0_DP + rim_distance_delta
         floor_distance=floor_diameter_old/2.0_DP+floor_distance_delta
         ejecta_distance=crater%rimwidth+rim_diameter_old/2.0_DP+ejecta_distance_delta
         rim_elevation=rim_elevation_delta+rim_height_old

         ! The following will be used in crater_profile:
         crater%fcrat=rim_distance*2.0_DP
         crater%floordiam=floor_distance*2.0_DP
         crater%rimheight=rim_elevation     
         
         if (radial_distance>ejecta_distance) then 
            newdem=0
            elchange  = newdem-surf(xpi,ypi)%dem       
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         endif

         if (radial_distance>rim_distance .and. radial_distance<ejecta_distance) then 
            newdem=crater_profile(user,crater,radial_distance/rim_distance)
            elchange  = newdem-surf(xpi,ypi)%dem       
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         endif

         if (radial_distance<rim_distance .and. radial_distance>floor_distance) then 
            newdem=crater_profile(user,crater,radial_distance/rim_distance)
            elchange  = newdem-surf(xpi,ypi)%dem       
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         endif

         if (radial_distance<floor_distance) then 
            newdem=-crater%floordepth
            elchange  = newdem-surf(xpi,ypi)%dem       
            deltaMtot = deltaMtot + elchange
            surf(xpi,ypi)%dem = newdem
         endif
      end do
   end do

   crater%rimheight =rim_height_old
   crater%floordiam=floor_diameter_old
   crater%fcrat=rim_diameter_old

   print *,crater%fcrat,crater%floordepth,crater%rimheight,crater%floordiam

end subroutine make_realistic_crater



 ! Universal shape model from PSD
subroutine Calculate_am_wl_phase_from_diameter(psd_1D,amplitude,wavelength,phase)
   use module_globals
   use module_realistic, EXCEPT_THIS_ONE => Calculate_am_wl_phase_from_diameter
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
   use module_realistic, EXCEPT_THIS_ONE => Calculate_breakpoint_slope_from_diameter
   use module_util
   implicit none
   ! in and out
   type(psdtype),intent(inout) :: psd_1D
   ! internal
   real(DP) :: random_number_normal
   
   ! excutable
   if (  psd_1D%diameter_in_km  <  psd_1D%diameter_in_km_trans  ) then
      psd_1D%bp2_x=psd_1D%bp2_x_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp2_x_k_b_sigma(2)         !+ random_number_normal * psd_1D%bp2_x_k_b_sigma(3) 
      psd_1D%bp2_y=psd_1D%bp2_y_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp2_y_k_b_sigma(2)     !+ random_number_normal * psd_1D%bp2_y_k_b_sigma(3) 
      psd_1D%bp3_y=psd_1D%bp3_y_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp3_y_k_b_sigma(2)    !+ random_number_normal * psd_1D%bp3_y_k_b_sigma(3) 
      psd_1D%bp4_y=psd_1D%bp4_y_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%bp4_y_k_b_sigma(2)   !+ random_number_normal * psd_1D%bp4_y_k_b_sigma(3)
      psd_1D%slope12=psd_1D%slope12_k_b_sigma(1)*psd_1D%diameter_in_km+psd_1D%slope12_k_b_sigma(2)   ! + random_number_normal * psd_1D%slope12_k_b_sigma(3)
      ! call util_random_number_normal(random_number_normal)

   else
      psd_1D%bp2_x=psd_1D%bp2_x_k_b_sigma(4)*psd_1D%diameter_in_km+psd_1D%bp2_x_k_b_sigma(5)!+ random_number_normal * psd_1D%bp2_x_k_b_sigma(6) 
      psd_1D%bp2_y=psd_1D%bp2_y_k_b_sigma(4)*psd_1D%diameter_in_km+psd_1D%bp2_y_k_b_sigma(5)!+ random_number_normal * psd_1D%bp2_y_k_b_sigma(6) 
      psd_1D%bp3_y=psd_1D%bp3_y_k_b_sigma(4)*psd_1D%diameter_in_km+psd_1D%bp3_y_k_b_sigma(5)!+ random_number_normal * psd_1D%bp3_y_k_b_sigma(6) 
      psd_1D%bp4_y=psd_1D%bp4_y_k_b_sigma(4)*psd_1D%diameter_in_km+psd_1D%bp4_y_k_b_sigma(5)!+ random_number_normal * psd_1D%bp4_y_k_b_sigma(6)
      psd_1D%slope12=psd_1D%slope12_k_b_sigma(4)*psd_1D%diameter_in_km+psd_1D%slope12_k_b_sigma(5)!+ random_number_normal * psd_1D%slope12_k_b_sigma(6)
      ! call util_random_number_normal(random_number_normal)
   endif
   !psd_1D%bp3_y=7
   print *,'bp3_y: ',psd_1D%bp3_y
   print *,'bp4_y: ',psd_1D%bp4_y
   call util_random_number_normal(random_number_normal)
   print *,'1',random_number_normal

end subroutine Calculate_breakpoint_slope_from_diameter

subroutine Calculate_targetPSD_from_breakpoint_slope(psd_1D,wavelength,psd)
   use module_globals
   use module_realistic, EXCEPT_THIS_ONE => Calculate_targetPSD_from_breakpoint_slope
   use module_util
   implicit none
   !in and out
   type(psdtype), intent(in)  :: psd_1D
   real(DP) ,dimension(:),allocatable,intent(out) :: wavelength,psd 
   ! internal
   integer(I4B) ::  i,bp2_x_index
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

      ! 3 segments, 4 breakpoints
   bp3_x=log10(psd_1D%input_x_max/2.0_DP) 
   bp4_x=log10(psd_1D%input_x_max) 

   do i = 1, psd_1D%num_sine
      bp2_x_index = i
      if (wavelength(i) > 10**psd_1D%bp2_x) exit
   end do
   !----------------------------------------------------------------------------------------------------------------------
   slope_12=psd_1D%slope12
   intercept_12 = psd_1D%bp2_y - slope_12 * psd_1D%bp2_x 

   slope_23=(psd_1D%bp3_y-psd_1D%bp2_y)/(bp3_x-psd_1D%bp2_x)
   intercept_23 = psd_1D%bp3_y - slope_23 * bp3_x 

   slope_34=(psd_1D%bp4_y-psd_1D%bp3_y)/(bp4_x-bp3_x)
   intercept_34 = psd_1D%bp4_y - slope_34 * bp4_x 
   !----------------------------------------------------------------------------------------------------------------------
   do i = 1, bp2_x_index-1
      call util_random_number_normal(random_number_normal)
      psd(i) = 10** (   slope_12*log10(wavelength(i))    +  intercept_12  )!+  random_number_normal*psd_1D%misfit_k_b_sigma  )      
   end do     
   do i = bp2_x_index, psd_1D%num_sine-2
      call util_random_number_normal(random_number_normal)
      psd(i) = 10** (   slope_23*log10(wavelength(i))    +  intercept_23 )! +  random_number_normal*psd_1D%misfit_k_b_sigma  )
   end do     
    
   call util_random_number_normal(random_number_normal)
   psd(psd_1D%num_sine-1) = 10** (   psd_1D%bp3_y )! +  random_number_normal*psd_1D%misfit_k_b_sigma  )
   call util_random_number_normal(random_number_normal)
   psd(psd_1D%num_sine) = 10** (   psd_1D%bp4_y )! +  random_number_normal*psd_1D%misfit_k_b_sigma  )
end subroutine Calculate_targetPSD_from_breakpoint_slope



subroutine Calculate_am_wl_phase_from_targetPSD(psd_1D,wavelength,psd,amplitude,phase)
   use module_globals
   use module_realistic, EXCEPT_THIS_ONE => Calculate_am_wl_phase_from_targetPSD
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


 

subroutine Make_outline(arc_length,psd_1D,amplitude,wavelength,phase,noise_component)
   use module_globals
   use module_realistic, EXCEPT_THIS_ONE => Make_outline
   implicit none
   ! in and out
   real(DP),intent(in) :: arc_length
   type(psdtype),intent(in) :: psd_1D
   real(DP),dimension(:),intent(in) :: amplitude
   real(DP),dimension(:),intent(in) :: wavelength
   real(DP),dimension(:),intent(in) :: phase
   real(DP),intent(out) :: noise_component
   ! internal
   integer(I4B) :: i
   real(DP) :: noise_component_ind


   ! print *,"aaaaaaaaa"


   ! print *,psd_1D%num_sine



   ! print *,0,amplitude(1) ,wavelength(1),phase(1)

   ! print *,1000,amplitude(1000) ,wavelength(1000),phase(1000)


   ! excutable
   noise_component = 0.0_DP
   do i = 1, psd_1D%num_sine

      ! print *,i,amplitude(i) ,wavelength(i),phase(i)
      noise_component_ind= amplitude(i) * sin( 2 * PI* 1/wavelength(i) * (  arc_length-   phase(i)   )     )      
      noise_component=noise_component+noise_component_ind
   end do
   noise_component=noise_component*1000.0_DP

end subroutine Make_outline
  
