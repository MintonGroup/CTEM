!****f* util/util_pop
! Name
!   util_pop -- Pop off a top layer. 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   
!   call util_pop(regolayer,oldregodata)
!
! DESCRIPTION
!    
!   Pop subroutine is to pop off the top head from an input old linked list.
!   This subroutine will be: 
!   * having a new pointer to point to the top of an old input layer. 
!   * if the second top is associated, linking the second top to be the new top.
!   * deallocating the pointer that inherits the original input layer. 
!
! ARGUMENTS
!   Input
!   * regolayer   -- pointer to the top of the regolith stack
!   
!   Output
!   * regolayer   -- pointer to the second top of the regolith stack
!   * oldregodata -- old regodta that was popped from the top of an input old layer
! 
! NOTES
!   If we reach the bottom of a layer, the layer that has a huge thickness is deallocated.
!   But this is considered a bug in the linked list related operation.
!
!***

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
   type(regolisttype),pointer :: oldhead

   ! Executable code

   if (associated(regolayer)) then
      oldhead  => regolayer 
      oldregodata = oldhead%regodata
      if (associated(oldhead%next)) then
         regolayer => oldhead%next
         deallocate(oldhead)
      else
         deallocate(regolayer)
         write(*,*) "util_pop error: We've reached the bottom of the regolith list!"
      end if
   else
      write(*,*) "util_pop error: Major error! The surface layer is not associated!"
   end if

   return
end subroutine util_pop

