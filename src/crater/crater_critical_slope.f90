!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Forms the maximum slope for the initial crater morphology
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


function crater_critical_slope(user,crater,iradsq) result(critical)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_form_exterior_func
   implicit none

  ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   integer(I4B),intent(in) :: iradsq
   real(DP) :: critical

   ! Internal variables

   critical = CRITSLP

end function crater_critical_slope

