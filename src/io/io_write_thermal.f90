!**********************************************************************************************************************************
!
!  Unit Name   : io_write_thermal
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for thermal grid
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine io_write_thermal(thermal)
    use module_globals
    use module_io, EXCEPT_THIS_ONE => io_write_thermal
    implicit none

    ! Arguments
    type(thermaltype),dimension(:,:,:),intent(in) :: thermal

    ! Internal variables
    integer(I4B), parameter :: LUN = 7

    !Executable code
    open(LUN,file=THERMFILE,status='replace',form='unformatted')

    write(LUN) thermal(:,:,:)%temperature

    close(LUN)

    return
end subroutine io_write_thermal


