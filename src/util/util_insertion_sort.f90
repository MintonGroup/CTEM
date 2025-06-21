!**********************************************************************************************************************************
!
!  Unit Name   : util_insertion_sort
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Array sort
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  Designed for time and temperature arrays of regolith layers
!
!**********************************************************************************************************************************
subroutine util_insertion_sort(arr)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_insertion_sort
    implicit none

    !Arguments
    real(SP), intent(inout) :: arr(:)

    !Internal variables
    integer :: i, j
    real(SP) :: key

    ! Executable code
    do i = 2, size(arr)
        key = arr(i)
        j = i - 1
        do while (j >= 1 .and. arr(j) > key)
            arr(j+1) = arr(j)
            j = j - 1
        end do
        arr(j+1) = key
    end do
end subroutine util_insertion_sort