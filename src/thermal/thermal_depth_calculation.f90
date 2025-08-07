!**********************************************************************************************************************************
!
!  Unit Name   : thermal_depth_calculation.f90
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculate depth and geotherm of each layer in the thermal grid relative to the DEM 
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       : The uppermost elevation must be the baseline for the array 
!
!**********************************************************************************************************************************
subroutine thermal_depth_calculation(user,surf,crater,domain,thermal,times,losses)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_depth_calculation
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(surftype),dimension(:,:),intent(in) :: surf
    type(cratertype),intent(in) :: crater
    type(domaintype),intent(inout) :: domain
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    real(DP),dimension(:,:,:),intent(inout) :: times,losses

    ! Internal variables
    integer(I4B) :: i,j,k, dd, sd, md, n
    real(DP) :: hmax, h, surfdepth, old_hmax, depth_difference, max_depth_difference
    type(thermaltype),dimension(:,:,:),allocatable :: old
    integer(I4B),dimension(:,:),allocatable :: shift
    real(DP),dimension(:,:,:),allocatable :: old_times, old_losses

    allocate(old_times,source=times)
    allocate(old_losses,source=losses)

    allocate(old,source=thermal(:,:,:))
    allocate(shift(user%gridsize,user%gridsize))
    do j=1,user%gridsize
        do i=1,user%gridsize
            do k=1,user%zgridsize
                if (thermal(i,j,k)%depth > 0) then
                    shift(i,j) = k
                    exit
                else
                    continue
                end if
            end do
        end do
    end do
    

    ! Executable Code

    hmax = maxval(surf(:,:)%dem)
    old_hmax = domain%hmax
    max_depth_difference = hmax - old_hmax
    md = max_depth_difference / user%zpix !this should be an integer pixel value

    do j=1,user%gridsize
        do i=1,user%gridsize
            h = surf(i,j)%dem
            n = 0
            do k=1,user%zgridsize
                surfdepth = hmax - h
                sd = int(surfdepth / user%zpix + 0.5)
                thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                thermal(i,j,k)%elevation = hmax - thermal(i,j,k)%relative_depth
                !if (thermal(i,j,k)%relative_depth < surfdepth) then !voxel is empty space above the surface
                if (k < sd) then
                    thermal(i,j,k)%temperature = 0.0_DP
                    thermal(i,j,k)%background = 0.0_DP
                    thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                    losses(i,j,k) = -1.0_DP
                    times(i,j,k) = -1.0_DP
                else
                    ! shift the temperature values by the difference between the new and old depth
                    depth_difference = thermal(i,j,k)%depth - old(i,j,k)%depth
                    if (depth_difference == 0.0_DP) exit
                    dd = nint(depth_difference / user%zpix)
                    if (dd == 0) exit
                    if (shift(i,j)+n .gt. user%zgridsize) then !temperature is equal to the background of the deepst voxel
                        thermal(i,j,k)%temperature = thermal(i,j,user%zgridsize)%background
                        losses(i,j,k) = -1.0_DP
                        times(i,j,k) = -1.0_DP
                    else if (shift(i,j)+n .lt. 1) then !? 
                        thermal(i,j,k)%temperature = 0.0_DP !adding new space (this should never be triggered because of the earlier check)
                        losses(i,j,k) = -1.0_DP
                        times(i,j,k) = -1.0_DP
                    else
                        thermal(i,j,k)%temperature = old(i,j,shift(i,j)+n)%temperature
                        losses(i,j,k) = old_losses(i,j,shift(i,j)+n)
                        times(i,j,k) = old_times(i,j,shift(i,j)+n)
                    end if
                    !thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                    thermal(i,j,k)%background = (13._DP/1000._DP) * thermal(i,j,k)%depth !13K/km for now; must match init_thermal.f90 (this should probably be a global variable)
                    n = n + 1
                end if             
            end do
        end do
    end do

    domain%hmax = hmax
    deallocate(old,shift,old_times,old_losses)

    return

end subroutine thermal_depth_calculation


