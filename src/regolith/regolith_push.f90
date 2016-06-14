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
subroutine regolith_push(surfi,newregodata,popflagi)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_push
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surfi
   type(regodatatype),intent(in) :: newregodata
   INTEGER(I4B),intent(inout) :: popflagi

   ! Internal variables
   type(regolisttype),POINTER :: newsurfi => null()
   integer(I4B):: allocstat

   ! Executable code
   !=======================================
   ! Initilize the linked list 
   !=======================================
   ! No value in the list 
   !=======================================
   IF (ASSOCIATED(surfi%regolayer)) THEN
      ALLOCATE(newsurfi, STAT=allocstat)
      IF (allocstat == 0) THEN
         NULLIFY(newsurfi%next)
         newsurfi%next => surfi%regolayer
         surfi%regolayer => newsurfi
         newsurfi%regodata = newregodata
         !IF (popflagi == 1) popflagi = 2
      ELSE
         PRINT *, 'Exhausted memory!'
      END IF
   END IF
   return
end subroutine regolith_push

