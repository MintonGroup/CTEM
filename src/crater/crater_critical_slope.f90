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


function crater_critical_slope(user,crater,lrad) result(critical)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_critical_slope
   implicit none

  ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: lrad
   real(DP) :: critical

   ! Internal variables
   real(DP) :: r

   r = lrad / crater%frad 
   if (r < 0.2_DP) then
      critical = 0.0_DP
   else if (r < 0.98_DP) then
      critical =  max(abs(0.228_DP + 2 * 0.083_DP * r - 3 * 0.039_DP * r**2),CRITSLP)
   else 
      critical =  max(abs(0.187_DP - 2 * 0.018_DP * r - 3 * 0.015_DP * r**2),CRITSLP)
   end if

end function crater_critical_slope

