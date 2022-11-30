!**********************************************************************************************************************************
!
!  Unit Name   : util_init_array
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initialize an new allocatable array
!  
!
!  Input
!    Arguments : regolayer :: array
!                oldlayer  :: old layer to pop off of the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       : 
!
!**********************************************************************************************************************************
subroutine util_init_array(regolayer,initstat)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_init_array
    implicit none
 
    ! Arguments
    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
    logical, intent(out)       :: initstat
 
    ! Internal variables
    integer(I4B) :: allocstat 
 
    ! Executable code
    initstat = .false.
    ! if (.not. associated(regolayer)) then
    !    allocate(regolayer, STAT=allocstat)
    !    if (allocstat == 0) then
    !       initstat = .true.
    !       nullify(regolayer%next)
          ! regolayer%regodata%thickness = sqrt(VBIG) ! This generates a buffer layer that the model should never reach if the run is structured properly
          ! regolayer%regodata%comp = 0.0_DP
          ! regolayer%regodata%meltfrac = 0.0_DP
          ! regolayer%regodata%porosity = 0.0_DP
          ! regolayer%regodata%age(:)   = 0.0_SP
    if (allocated(regolayer)) deallocate(regolayer)
    allocate(regolayer(1))
    regolayer(1)%thickness = sqrt(VBIG) ! This generates a buffer layer that the model should never reach if the run is structured properly
    regolayer(1)%comp = 0.0_DP
    regolayer(1)%meltfrac = 0.0_DP
    regolayer(1)%porosity = 0.0_DP
    regolayer(1)%age(:)   = 0.0_SP
    !    else
    !       write(*,*) 'util_init_list: Initialization failed. Exhausted memory.'
    !    end if
    ! else
    !    write(*,*) 'util_init_list: Initialization failed. Regolayer already associated.'
    ! end if
 
    return
 end subroutine util_init_array