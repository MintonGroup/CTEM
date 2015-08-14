!**********************************************************************************************************************************
!
!  Unit Name   : regolith_push
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pushes a new regolith block onto an old surface
!  
!
!  Input
!    Arguments : regolayer :: pointer to the top of the regolith stack
!                newlayer  :: new layer to push onto the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_push(toplayer,newlayer)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_push
   implicit none

   ! Arguments
   type(regolayertype),pointer :: toplayer
   type(regolayertype),intent(in) :: newlayer

   ! Internal variables
   type(regolayertype),pointer :: current

   ! Executable code
   current = toplayer
   current%next => toplayer
   toplayer => current

   return
end subroutine regolith_push

