!****f* util/util_init_list
! Name
!   util_init_list -- Initialize an new linked list (see DESCRIPTION).
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   
!   call util_init_list(regolayer,initstat)
!
! DESCRIPTION
!    
!   This initialization process to a new linked list is to
!   * Check the associated state of a new linked list. Idealiy, it should not be associated. 
!   * Nullify the head of this new linked list.
!   * Allocate a space (container) for regodata-typed data to the head.
!   * Assign corresponding values to each regodata type.
!
! ARGUMENTS
!   Input
!   * regolayer -- pointer to the top of the regolith stack
!   
!   Output
!   * initstat  -- return logical, T or F, about the success of initializing a list. 
! 
! NOTES
!   The initial thickness is a very huge number that makes sure the "pop" never gets to it.
!   This subroutine suggests that an initialized "popped" list has the huge number at the last layer, 
!   so please be cautious about any calculations involving popped list or using olddata that is popped
!   and gets collected by popped list should avoid using the value of this last layer. 
!
!***

!**********************************************************************************************************************************
!
!  Unit Name   : util_init_list
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initialize an new list
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
!  Notes       : This subroutine suggests that a popped list has this value at the last layer. Any calculation involving popped list
!                or using olddata that is popped and gets collected by popped list should avoid using the value of the last
!                layer. 
!
!**********************************************************************************************************************************
subroutine util_init_list(regolayer,initstat)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_init_list
   implicit none

   ! Arguments
   type(regolisttype),pointer :: regolayer
   logical, intent(out)       :: initstat

   ! Internal variables
   integer(I4B) :: allocstat 

   ! Executable code
   initstat = .false.
   if (.not. associated(regolayer)) then
      allocate(regolayer, STAT=allocstat)
      if (allocstat == 0) then
         initstat = .true.
         nullify(regolayer%next)
         regolayer%regodata%thickness = sqrt(VBIG) ! This generates a buffer layer that the model should never reach if the run is structured properly
         regolayer%regodata%comp = 0.0_DP
         regolayer%regodata%meltfrac = 0.0_DP
         regolayer%regodata%porosity = 0.0_DP
         regolayer%regodata%age(:)   = 0.0_SP
      else
         write(*,*) 'util_init_list: Initialization failed. Exhausted memory.'
      end if
   else
      write(*,*) 'util_init_list: Initialization failed. Regolayer already associated.'
   end if

   return
end subroutine util_init_list

