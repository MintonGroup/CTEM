!**********************************************************************************************************************************
!
!  Unit Name   : regolith_pop
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pushes a old regolith block onto an old surface
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
subroutine regolith_pop(surf)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_pop
   implicit none

   ! Arguments
   type(surftype),intent(inout):: surf

   ! Internal variables
   type(regolayertype),pointer :: current

   ! Executable code
!   if (.not. associated(surf%regolayer%next)) then
!      write(*,*) "Error: Dug too deep. Beware of balrog."
!   else
!      deallocate (surf%regolayer%next)
!      surf%regolayer => surf%regolayer%next
!   end if

   current => surf%regolayer
   if (.not. associated(surf%regolayer%next)) then
      write(*,*) "Error: Dug too deep. Beware of balrog."
   else
      surf%regolayer => surf%regolayer%next
      !deallocate (current)
   end if
   return
end subroutine regolith_pop

