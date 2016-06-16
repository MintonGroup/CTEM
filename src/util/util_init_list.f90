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
!  Notes       :  
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
      else
         write(*,*) 'util_init_list: exhausted memory.'
      end if
   else
      write(*,*) 'util_init_list: Initialization went wrong. regolayer already associated.'
   end if
   regolayer%regodata%thickness = VBIG ! This generates a buffer layer that the model should never reach if the run is structured properly
   regolayer%regodata%comp = 0.0_DP

   return
end subroutine util_init_list

