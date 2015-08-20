!**********************************************************************************************************************************
!
!  Unit Name   : regolith_rays
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Determine a history of a stream tube
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_rays(user,crater,domain,ejtble,ejb)
   use module_globals 
   use module_ejecta
   use module_regolith, EXCEPT_THIS_ONE => regolith_rays
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   integer(I4B),intent(in) :: ejtble
   type(ejbtype),dimension(ejtble),intent(in)   :: ejb

   ! Traversing a linked list 
   real(DP),parameter :: a = 0.936457 ! Fitting parameters for the relation between height difference and a radial position of a stream tube
   real(DP),parameter :: b = 1.12368
   real(DP),parameter :: deltar = 0.1_DP
   real(DP)     :: outeredge,inneredge,rbead,ce1sq,excav,landarea
   real(DP)     :: vseg,theta,thick,area_ray,area_tot
   integer(I4B) :: i,j,nloop,nbead
   real(DP)     :: lrad1,vejsq1,ejang1
   real(DP)     :: lrad2,vejsq2,ejang2
   logical      :: firstrun
   
   ! Executalbe code

   ! The range of crater ray bead modeling
   outeredge = ejb(1000)%erad !crater%frad 
   inneredge = ejb(1)%erad
   nloop = int(abs(outeredge - inneredge)/ (2.0 * deltar))
   area_ray = 0._DP
   area_tot = 0._DP
   write(*,*) outeredge,inneredge,nloop 
   do i=1,nloop-1
      ! Location of simulated ray bead within a final crater 
      rbead = outeredge + 2.0 * deltar * (real(i)-0.5)
      ! How many simulated ray bead is estimated, and will be used to estimate the area of a simulated ray at the landing site. 
      ! Compare to the ejecta blanket thickness at a landing site, the area that a simulated ray bead can be obtained. 
      nbead = int(2.0 * PI * rbead / (2.0 * deltar))
      firstrun = .true.
      do j=1,1!nbead
         ! The simulated ray bead's direction 
         theta = 2.0 * deltar /rbead * (real(j)-0.5) / PI * 180.0
         ! The total volume of a simulated ray bead
         vseg = 0.25_DP * PI * deltar**2 * a**2 * rbead / b * abs(tan(b) - b) + sqrt(2.0_DP)/2.0_DP*PI*deltar**3
         ! Calculate the total volume of ejecta blanket thickness in CTEM
         call ejecta_blanket(user,crater,domain,rbead-deltar,lrad1,vejsq1,ejang1,firstrun)
         firstrun = .true.   
         call ejecta_blanket(user,crater,domain,rbead+deltar,lrad2,vejsq2,ejang2,firstrun)
         !call ejecta_thickness(user,crater,rbead-deltar,rbead+deltar,lrad1,lrad2,thick)
         !area_tot = PI * abs(lrad1**2 - lrad2**2)
         ce1sq = (2._DP / (CT * CT)) * user%mu_r**2 / (1._DP + 2 * user%mu_r + user%mu_r**2)
         excav = abs(0.5_DP * ce1sq * PI * abs((rbead+deltar)**3 - (rbead-deltar)**3))
         landarea = PI * abs(lrad2*lrad2 - lrad1*lrad1)
         thick = excav / landarea
         area_ray = vseg / thick 
      end do
      write(*,*) (lrad1+lrad2)/2.0,area_ray/(2.0*deltar)/(PI * (lrad1+lrad2)) * 180.0 
   end do
  
   return
end subroutine regolith_rays
