!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the degradation function from Minton et al. (2019)
!  
!
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


function crater_degradation_function(user,r) result(Kd)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_degradation_function
   implicit none

  ! Arguments
   type(usertype),intent(in) :: user
   real(DP),intent(in) :: r
   real(DP) :: Kd
   real(DP) :: A,r_break,alpha_1,alpha_2,delta

   ! Testing values that match both mare scale and highlands scale
   alpha_1 = user%psi
   alpha_2 = user%psi2

   r_break = user%rbreak
   delta = 1.0_DP

   A = user%Kd1 * r_break**(alpha_1) / (1_DP + (1.d0 / delta) * (alpha_1 - alpha_2) / alpha_1)**2 
   Kd =  A * (r / r_break)**(alpha_1) * (0.5_DP * (1.0_DP + (r / r_break)**(1.0_DP / delta)))**((alpha_2 - alpha_1 ) / delta)

   return

end function crater_degradation_function

