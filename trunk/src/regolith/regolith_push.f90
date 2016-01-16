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
subroutine regolith_push(surf,newlayer)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_push
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surf
   type(regolayertype),intent(in) :: newlayer

   ! Internal variables
   type(regolayertype),pointer :: current
   integer(I4B):: err

   ! Executable code
   ! make sure if available memory for a new node
   allocate(current, stat = err) 

!   if (newlayer%thickness = 0._DP) then 
!      surf%regolayer%thickness = surf%regolayer%thickness
!   end if

   if (err == 0) then
      nullify(current%next)                     ! initialize the pointer of a new node
      current%thickness = newlayer%thickness
      current%meltfrac  = newlayer%meltfrac
      current%comp = newlayer%comp
      current%next   => surf%regolayer
      surf%regolayer => current
   else
      write(*,*) 'exhausted memory.'
   end if

   return
end subroutine regolith_push

