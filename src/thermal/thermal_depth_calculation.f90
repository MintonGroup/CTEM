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
subroutine thermal_depth_calculation(user,surf,crater,domain,thermal)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_depth_calculation
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(surftype),dimension(:,:),intent(in) :: surf
    type(cratertype),intent(in) :: crater
    type(domaintype),intent(inout) :: domain
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal

    ! Internal variables
    integer(I4B) :: i,j,k, dd, sd, md
    real(DP) :: hmax, h, surfdepth, old_hmax, depth_difference, max_depth_difference
    type(thermaltype),dimension(:,:,:),allocatable :: old

    allocate(old,source=thermal(:,:,:))
    

    ! Executable Code

    hmax = maxval(surf(:,:)%dem)
    old_hmax = domain%hmax
    max_depth_difference = hmax - old_hmax
    md = max_depth_difference / user%zpix !this should be an integer pixel value

    do j=1,user%gridsize
        do i=1,user%gridsize
            h = surf(i,j)%dem
            do k=1,user%zgridsize
                surfdepth = hmax - h
                sd = int(surfdepth / user%zpix + 0.5)
                thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - hmax
                thermal(i,j,k)%elevation = hmax - thermal(i,j,k)%relative_depth
                !if (thermal(i,j,k)%relative_depth < surfdepth) then !voxel is empty space above the surface
                if (k < sd) then
                    thermal(i,j,k)%temperature = 0.0_DP
                    thermal(i,j,k)%background = 0.0_DP
                    thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                else
                    ! shift the temperature values by the difference between the new and old depth
                    depth_difference = thermal(i,j,k)%depth - old(i,j,k)%depth
                    dd = nint(depth_difference / user%zpix)
                    if (k-dd .gt. user%zgridsize) then !temperature is equal to the background of the deepst voxel
                        thermal(i,j,k)%temperature = thermal(i,j,user%zgridsize)%background
                    else if (k-dd .lt. 1) then !? 
                        thermal(i,j,k)%temperature = 0.0_DP !adding new space (this should never be triggered because of the earlier check)
                    else
                        thermal(i,j,k)%temperature = old(i,j,k-dd)%temperature
                    end if
                    !thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                    thermal(i,j,k)%background = (13._DP/1000._DP) * thermal(i,j,k)%depth !13K/km for now; must match init_thermal.f90 (this should probably be a global variable)
                end if             
            end do
        end do
    end do

    domain%hmax = hmax
    deallocate(old)

    return

end subroutine thermal_depth_calculation


