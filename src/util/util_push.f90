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
subroutine util_push(surfi,newregodata,popflagi)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_push
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surfi
   type(regodatatype),intent(in) :: newregodata
   integer(I4B),intent(inout) :: popflagi

   ! Internal variables
   type(regolisttype),POINTER :: newsurfi => null()
   integer(I4B):: allocstat

   ! Executable code
   !=======================================
   ! Initilize the linked list 
   !=======================================
   ! No value in the list 
   !=======================================
   if (associated(surfi%regolayer)) then
      allocate(newsurfi, stat=allocstat)
      if (allocstat == 0) then
         nullify(newsurfi%next)
         newsurfi%next => surfi%regolayer
         surfi%regolayer => newsurfi
         newsurfi%regodata = newregodata
         !IF (popflagi == 1) popflagi = 2
      else
         write(*,*) 'Exhausted memory!'
      end if
   end if
   return
end subroutine util_push

