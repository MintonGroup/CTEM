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
subroutine util_push_regotemp(regotemp,regotime)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_push_regotemp
    implicit none

    ! Arguments
    real(SP),dimension(:,:),allocatable,intent(inout) :: regotemp
    real(SP),dimension(:,:),allocatable,intent(inout) :: regotime


    ! Internal variables
    real(SP), dimension(:,:), allocatable :: newtemp
    real(SP), dimension(:,:), allocatable :: newtime

    integer(I4B) :: nold, NC

    ! Executable code

    if (.not. allocated(regotemp)) then
       allocate(regotemp(1,1))
       regotemp(1,1) = 0.0_SP
    end if

    if (.not. allocated(regotime)) then
        allocate(regotime(1,1))
        regotime(1,1) = 0.0_SP
    end if

    NC = size(regotemp(:,1))
    nold = size(regotemp(1,:))

    allocate(newtemp(NC,nold+1))
    newtemp(:,1:nold) = regotemp(:,1:nold)
    newtemp(:,nold+1) = 0.0_SP ! This will be filled in the thermal_link subroutine
    call move_alloc(newtemp,regotemp)


    allocate(newtime(NC,nold+1))
    newtime(:,1:nold) = regotime(:,1:nold)
    newtime(:,nold+1) = 0.0_SP ! This will be filled in the thermal_link subroutine
    call move_alloc(newtime, regotime)
end subroutine util_push_regotemp