!**********************************************************************************************************************************
!
!  Unit Name   : util_push_regotemp
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pushes a new regolith temperarture structure
!  
!
!  Input
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine util_push_regotemp(regolayer)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_push_regotemp
    implicit none

    ! Arguments
    type(regodatatype),intent(inout) :: regolayer

    ! Internal variables
    real(SP), dimension(:), allocatable :: newtemp
    real(SP), dimension(:), allocatable :: newtime

    integer(I4B) :: nold, NC

    ! Executable code
    if (.not. allocated(regolayer%regotemp)) then
       allocate(regolayer%regotemp(1))
       regolayer%regotemp(1) = 0.0_SP
    end if

    if (.not. allocated(regolayer%regotime)) then
        allocate(regolayer%regotime(1))
        regolayer%regotime(1) = 0.0_SP
    end if

    nold = size(regolayer%regotemp)

    allocate(newtemp(nold+1))
    newtemp(1:nold) = regolayer%regotemp(1:nold)
    newtemp(nold+1) = 0.0_SP ! This will be filled in the thermal_loss_link subroutine
    call move_alloc(newtemp,regolayer%regotemp)


    allocate(newtime(nold+1))
    newtime(1:nold) = regolayer%regotime(1:nold)
    newtime(nold+1) = 0.0_SP ! This will be filled in the thermal_loss_link subroutine
    call move_alloc(newtime,regolayer%regotime)
end subroutine util_push_regotemp