!**********************************************************************************************************************************
!
!  Unit Name   : util_destroy_list
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Destroyes an old list
!  
!
!  Input
!    Arguments : regolayer :: pointer to the top of the regolith stack
!                oldlayer  :: old layer to pop off of the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine util_destroy_list(regolayer)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_destroy_list
   implicit none

   ! Arguments
   type(regolisttype),pointer :: regolayer

   ! Internal variables
   type(regodatatype) :: oldregodata

   ! Executable code

   do while(associated(regolayer))
      call util_pop(regolayer,oldregodata)
   end do

   return
end subroutine util_destroy_list

