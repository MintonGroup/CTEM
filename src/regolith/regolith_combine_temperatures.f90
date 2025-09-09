!**********************************************************************************************************************************
!
!  Unit Name   : regolith_combine_temperatures
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds maximum value of Ar loss in a single mixed layer      
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
! 
!  Notes       :  This has been "retrofitted" from a routine used to combine temperatures. As such, some variable names won't make sense (e.g. "regotemp" is actually Ar loss)
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
    real(SP), dimension(:), allocatable :: temp_times(:), max_temps(:)
    logical :: found
    real(SP), allocatable :: sum_temp(:), sum_weight(:), avg_temps(:)
    real(SP) :: t, weight, local_sum, local_cnt

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

! Step 3: Accumulate ALL matching entries (no exit), counting duplicates too
    do i = 1, N
        NT = size(layers(i)%regotime)
        do j = 1, M
            t = newlayer%regotime(j)
            local_sum = 0.0_SP
            local_cnt = 0.0_SP
    
            do k = 1, NT
                if (abs(layers(i)%regotime(k) - t) < 1.0e-6_SP) then
                    local_sum = local_sum + layers(i)%regotemp(k)
                    local_cnt = local_cnt + 1.0_SP
                end if
            end do
    
            ! Add ALL matches from this layer for time t
            if (local_cnt > 0.0_SP) then
                sum_temp(j)   = sum_temp(j)   + local_sum
                sum_weight(j) = sum_weight(j) + local_cnt
            end if
        end do
    end do
    
    ! Step 4: Compute averages
    allocate(avg_temps(M))
    do j = 1, M
        if (sum_weight(j) > 0.0_SP) then
            avg_temps(j) = sum_temp(j) / sum_weight(j)
        else
            avg_temps(j) = 0.0_SP 
        end if
    end do

    newlayer%regotemp = max_temps


    deallocate(sum_temp,sum_weight,temp_times,max_temps)

    return

end subroutine regolith_combine_temperatures