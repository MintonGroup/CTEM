!****f* util/util_push
! Name
!   util_push -- Push a new layer onto the top of an old layer. 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   
!   call util_push(regolayer,newregodata)
!
! DESCRIPTION
!    
!   Push subroutine is to create new head and push it with new data onto the old layer.
!   This subroutine will be: 
!   * checking if the head of an input old layer is associated.
!   * if we have enough space, allocating a new head and space for new data.  
!   * at the end, linking new head with new data to the old head. 
!
! ARGUMENTS
!   Input
!   * regolayer   -- pointer to the top of the regolith stack
!   * newregodata -- new regodata that is about to be pushed. 
!   
!   Output
!   * regolayer   -- pointer to the top of the regolith stack with newlayer. 
! 
! NOTES
! 
!***

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

