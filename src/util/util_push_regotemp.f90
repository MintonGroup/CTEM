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
subroutine util_push_regotemp(regotemp)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_push_regotemp
    implicit none

    ! Arguments
    type(thermalhisttype),dimension(:),allocatable,intent(inout) :: regotemp

    ! Internal variables
    type(thermalhisttype), dimension(:), allocatable :: newlayer
    integer(I4B) :: nold

    ! Executable code

    nold = size(regotemp)

    if (.not. allocated(regotemp)) then
        write(*,*) "ERROR"
    end if

    allocate(newlayer(nold+1))
    newlayer(1:nold) = regotemp(1:nold)
    newlayer(nold+1)%temperature = 0.0_DP ! This will be filled in the thermal_link subroutine
    newlayer(nold+1)%time = 0.0_DP ! This will be filled in the thermal_link subroutine
    newlayer(nold+1)%timeGa = 0.0_DP ! This will be filled in the thermal_link subroutine
    call move_alloc(newlayer, regotemp)
end subroutine util_push_regotemp