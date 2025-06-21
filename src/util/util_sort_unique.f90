!**********************************************************************************************************************************
!
!  Unit Name   : util_sort_unique
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Combine and sort arrays
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
subroutine util_sort_unique(input,output)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_sort_unique
    implicit none
    
    !Arguments
    real(SP),intent(in) :: input(:)
    real(SP),allocatable,intent(out) :: output(:)

    !Internal variables
    real(SP), allocatable :: sorted(:)
    integer(I4B) :: i, n, count

    sorted = input
    call util_insertion_sort(sorted)

    ! Remove duplicates
    n = size(sorted)
    allocate(output(n))
    count = 1
    output(count) = sorted(1)

    do i = 2, n
        if (abs(sorted(i) - sorted(i-1)) > 1.0e-6_SP) then
            count = count + 1
            output(count) = sorted(i)
        end if
    end do

    call move_alloc(output, sorted)
    allocate(output(count))
    output = sorted(:count)
    return
end subroutine util_sort_unique

