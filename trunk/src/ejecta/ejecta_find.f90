!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_find
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds the ejecta thickness at a given radius
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
subroutine ejecta_find(user,crater,domain,lrad,thick,strflag,ce1sq,ce2sq,bit,yterm,x1)
   use module_globals
   use module_util
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_find
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: ce1sq,ce2sq,bit,yterm,lrad
   real(DP),intent(out) :: thick
   integer(I4B),intent(in) :: strflag
   real(DP),intent(inout) :: x1

   ! Internals

   ! Executable code

   return
end subroutine ejecta_find

