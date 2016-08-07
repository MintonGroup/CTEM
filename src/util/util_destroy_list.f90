!****f* util/util_destroy_list
! Name
!   util_destroy_list -- Destroy a linked list (see DESCRIPTION).
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   
!   call util_destroy_list(regolayer)
!
! DESCRIPTION
!    
!   This subroutine takes an old linked list (layer) and poppes off all available tops that are stored in
!   this old layer. At the end, it releases the memory of this old linked list. 
!
! ARGUMENTS
!   Input
!   * regolayer -- pointer to the top of the regolith stack
!   
! NOTES
!   The association state of an output regolayer is not associated.
! 
!***

!**********************************************************************************************************************************
!
!  Unit Name   : util_destroy_list
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Destroyes an old list
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
subroutine util_destroy_list(regolayer)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_destroy_list
   implicit none

   ! Arguments
   type(regolisttype),pointer :: regolayer

   ! Internal variables
   type(regodatatype) :: oldregodata
   ! Executable code

   do while (associated(regolayer%next))
      call util_pop(regolayer,oldregodata)
   end do

   if (associated(regolayer)) then
      deallocate(regolayer)
   else
      write(*,*) 'util_destroy_list: the head not associated.'
   end if

   return
end subroutine util_destroy_list

