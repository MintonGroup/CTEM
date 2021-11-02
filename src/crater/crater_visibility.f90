!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the visibility function from Minton et al. (2019)
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


function crater_visibility(user,crater,Kval) result(iscountable)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_visibility
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: Kval

   ! Result variable
   logical :: iscountable 
   
   ! Internal variables
   real(DP) :: Kv
   real(DP),parameter :: Kv1 = 0.17
   real(DP),parameter :: gam = 2.0

   ! Counting parameter for simple craters from Bryan Howl's study in Minton et al. (2019)
   
   Kv = Kv1 * (crater%frad)**gam
   if (Kval >= kv) then
      iscountable = .false.
   else
      iscountable = .true.
   end if

   return

end function crater_visibility

