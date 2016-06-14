!**********************************************************************************************************************************
!
!  Unit Name   : util_pop
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
subroutine util_pop(surfi,popflagi)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_pop
   implicit none

   ! Arguments
   type(surftype),intent(inout):: surfi
   integer(I4B),intent(inout)  :: popflagi

   ! Internal variables
   type(regolisttype),pointer :: current  => null()

   ! Executable code
   !===========================================
   ! Check if the head is associated
   !===========================================
   !IF (.NOT. ASSOCIATED(surfi%regolayer)) RETURN
   !current  => surfi%regolayer
   !previous => current%next 
   !IF (ASSOCIATED(previous%next)) THEN 
   !PRINT *, current%regodata%thickness, previous%next%regodata%thickness, surfi%regolayer%regodata%thickness
   !previous => previous%next
   !surfi%regolayer => previous
   !PRINT *, previous%regodata%thickness, surfi%regolayer%regodata%thickness
   !DEALLOCATE(current)
   !END IF

   if (associated(surfi%regolayer)) then
      current  => surfi%regolayer 
      if (associated(current%next)) then
         surfi%regolayer => surfi%regolayer%next
         deallocate(current)
      end if
   end if

   return
end subroutine util_pop

