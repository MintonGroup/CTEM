!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_blanket_func
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Wrapper function for root finding the ejecta_blanket_func
!  subroutine by solving for erad as a function of lrad
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
function ejecta_blanket_func(user,crater,domain,erad,lrad,vejsq,ejang,firstrun) result(ans)
   use module_globals
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_blanket_func
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: erad,lrad
   logical,intent(inout) :: firstrun
   real(DP),intent(out) :: vejsq,ejang
   real(DP) :: ans

   ! Internals
   real(DP) :: lradguess

   ! Executable code
   lradguess = lrad
   call ejecta_blanket(user,crater,domain,erad,lradguess,vejsq,ejang,firstrun)
   ans = lradguess - lrad

   return
end function ejecta_blanket_func

