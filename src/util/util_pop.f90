!**********************************************************************************************************************************
!
!  Unit Name   : util_pop
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Popped a old regolith block onto an old surface
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
subroutine util_pop(regolayer,oldregodata)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_pop
   implicit none

   ! Arguments
   type(regolisttype),pointer :: regolayer
   type(regodatatype),intent(out) :: oldregodata

   ! Internal variables
   type(regolisttype),pointer :: current  => null()

   ! Executable code

   if (associated(regolayer)) then
      current  => regolayer 
      oldregodata = current%regodata
      if (associated(current%next)) then
         regolayer => regolayer%next
         deallocate(current)
      else
         regolayer => null()
         write(*,*) "util_pop error: We've reached the bottom of the regolith list!"
      end if
   else
      write(*,*) "util_pop error: Major error! The surface layer is not associated!"
   end if

   return
end subroutine util_pop

