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
         regolayer%regodata%thickness = VBIG ! This generates a buffer layer that the model should never reach if the run is structured properly
         regolayer%regodata%comp = 0.0_DP
         regolayer%regodata%meltfrac = 0.0_DP
         regolayer%regodata%porosity = 0.0_DP
      else
         write(*,*) 'util_init_list: Initialization failed. Exhausted memory.'
      end if
   else
      write(*,*) 'util_init_list: Initialization failed. Regolayer already associated.'
   end if

   return
end subroutine util_init_list

