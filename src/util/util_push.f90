!**********************************************************************************************************************************
!
!  Unit Name   : util_push
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
subroutine util_push(regolayer,newregodata)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_push
   implicit none

   ! Arguments
   type(regolisttype),pointer :: regolayer
   type(regodatatype),intent(in) :: newregodata

   ! Internal variables
   type(regolisttype),pointer :: newlayer => null()
   integer(I4B):: allocstat

   ! Executable code
   !=======================================
   ! Initilize the linked list 
   !=======================================
   ! No value in the list 
   !=======================================
   if (associated(regolayer)) then
      allocate(newlayer, stat=allocstat)
      if (allocstat == 0) then
         nullify(newlayer%next)
         newlayer%next => regolayer
         newlayer%regodata = newregodata
         regolayer => newlayer
      else
         write(*,*) 'util_push error: Exhausted memory!'
      end if
   else
      write(*,*) "util_push error: regolayer is not associated!"
   end if
   return
end subroutine util_push

