!**********************************************************************************************************************************
!
!  Unit Name   : regolith_pop
!  Unit Type   : subroutine
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
subroutine regolith_pop(regolayer,oldlayer)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_pop
   implicit none

   ! Arguments
   type(regolayertype),pointer :: regolayer
   type(regolayertype),intent(out) :: oldlayer

   ! Internal variables

   ! Executable code
   if (.not.associated(regolayer%next)) then
      write(*,*) "Error: Dug too deep. Beware of balrog."
   else
      oldlayer = regolayer
      regolayer => regolayer%next
   end if


   return
end subroutine regolith_pop

