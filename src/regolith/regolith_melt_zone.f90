!**********************************************************************************************************************************
!
!  Unit Name   : regolith_melt_zone
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Computes the radius of melt zone
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_melt_zone(user,crater,dimp,vimp,rmelt,depthb)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_melt_zone
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   real(DP),intent(in)  :: dimp,vimp    ! diameter and impact velocity of projectile for testing case and a real run
   real(DP),intent(out) :: rmelt, depthb

   ! Internal variables
   real(DP),parameter      :: Em = 3.42d06 ! specific internal energy for highland (Bjorkman and Holsapple 1987)
   real(DP)                :: rimp
   real(DP)                :: volm,vtc
   real(DP)                :: b,c,d,e

   ! Executable code 

   rimp = 0.5 * dimp
   depthb = 0.5 * dimp
   vtc = PI/3.0 * (crater%rad)**3
   volm = 2.9 * vtc * Em**(-0.85) * (dimp)**0.66 * (user%gaccel) **0.66 * vimp**(0.37) 
   ! Calculate radius of melt zone (shifted melt zone model Pierazzo et al. 1997)
   ! rmelt = (rvapor**3 + 1.5/PI*volm)**(1.0/3.0) (hemisphere model)
   ! R_m ^3 + 1.5 * R_imp * R_m ^2 - 2.5 * R_imp ^3 - 3/(2 * PI) * volm = 0
   ! ax^3 + bx^2 + c + d = 0, where c = - 2.5 * R_imp ^3, d = - 3/(2 * PI) * volm
   b = 1.5 * depthb
   c = -0.5 * depthb**3 -2.0 * rimp**(3)
   d = -1.5/PI * volm
   e = ( sqrt(-4.0 * b **6 + (-2.0 * b**3 - 27.0 * c - 27.0 * d)**2) - 2.0 * b**3 - 27.0 * c - 27.0 * d)**(1.0/3.0)
   rmelt =  -1.0/3.0 * b + (2)**(1.0/3.0)*b**2/(3.0 * e) + e/(3.0 * (2)**(1.0/3.0))

   return
end subroutine regolith_melt_zone
