!**********************************************************************************************************************************
!
!  Unit Name   : regolith_combine_temperatures
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Combines temperatures in the layer structure (associated with the thermal part of mixing)        
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
subroutine regolith_combine_temperatures(newlayer,layers,N)
    use module_globals
    use module_util
    use module_regolith, EXCEPT_THIS_ONE => regolith_mix
    implicit none

    ! Arguments
    type(regodatatype),intent(inout) :: newlayer
    type(regodatatype),dimension(:),allocatable,intent(in) :: layers
    integer(I4B),intent(in) :: N

    ! Internal variables
    integer(I4B) :: i, j, k, total_times, count, NT, M
    real(SP), dimension(:), allocatable :: temp_times(:)
    logical :: found
    real(SP), allocatable :: sum_temp(:), sum_weight(:)
    real(SP) :: t, weight

    ! Executable code

    ! Step 1: Collect and deduplicate all time points
    total_times = 0
    do i = 1, N
        total_times = total_times + size(layers(i)%regotime)
    end do

    allocate(temp_times(total_times))
    count = 0
    do i = 1, N
        NT = size(layers(i)%regotime)
        do j = 1,NT
            count = count + 1
            temp_times(count) = layers(i)%regotime(j)
        end do
    end do


    call util_sort_unique(temp_times, newlayer%regotime)
    M = size(newlayer%regotime)

    ! Step 2: Initialize accumulators
    allocate(sum_temp(M), sum_weight(M))
    sum_temp = 0.0_SP
    sum_weight = 0.0_SP

    ! Step 3: Populate accumulators
    do i = 1, N
        weight = layers(i)%thickness
        do j = 1, M
            t = newlayer%regotime(j)
            found = .false.
            NT = size(layers(i)%regotemp)
            do k = 1, NT
                if (abs(layers(i)%regotime(k) - t) < 1.0e-6_SP) then
                    sum_temp(j) = sum_temp(j) + weight * layers(i)%regotemp(k)
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                sum_temp(j) = sum_temp(j) + weight * 0.0_SP !Background won't work as it didn't link
            end if
            sum_weight(j) = sum_weight(j) + weight
        end do
    end do

    ! Step 4: Compute weighted average
    allocate(newlayer%regotemp(M))
    do j = 1, M
        if (sum_weight(j) > 0.0_SP) then
            newlayer%regotemp(j) = sum_temp(j) / sum_weight(j)
        else
            newlayer%regotemp(j) = -1.0_SP  ! or mark as NaN
        end if
    end do

    deallocate(sum_temp,sum_weight,temp_times)

    return

end subroutine regolith_combine_temperatures