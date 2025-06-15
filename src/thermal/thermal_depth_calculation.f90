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
    integer(I4B) :: i,j,k, dd
    real(DP) :: hmax, h, surfdepth, old_hmax, depth_difference
    real(DP),dimension(:,:,:),allocatable :: oldtemps

    allocate(oldtemps,source=thermal(:,:,:)%temperature)
    

    ! Executable Code

    hmax = maxval(surf(:,:)%dem)
    old_hmax = domain%hmax
    depth_difference = hmax - old_hmax
    dd = depth_difference / user%zpix !this should be an integer pixel value

    do j=1,user%gridsize
        do i=1,user%gridsize
            h = surf(i,j)%dem
            do k=1,user%zgridsize
                surfdepth = hmax - h
                dd = surfdepth / user%zpix !This only works if the old surfdepth is 0
                thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - hmax
                thermal(i,j,k)%elevation = hmax - thermal(i,j,k)%relative_depth
                if (thermal(i,j,k)%relative_depth < surfdepth) then !voxel is empty space above the surface
                    thermal(i,j,k)%temperature = 0.0_DP
                    thermal(i,j,k)%background = 0.0_DP
                    thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                else
                    ! shift the temperature values by the difference between the new and old hmax
                    if (k-dd .gt. user%zgridsize) then !temperature is equal to the background of the deepst voxel
                        thermal(i,j,k)%temperature = thermal(i,j,user%zgridsize)%background
                    else if (k-dd .lt. 1) then !? 
                        thermal(i,j,k)%temperature = 0.0_DP !adding new space (this should never be triggered because of the earlier check)
                    else
                        thermal(i,j,k)%temperature = oldtemps(i,j,k-dd)
                    end if
                    thermal(i,j,k)%depth = thermal(i,j,k)%relative_depth - surfdepth
                    thermal(i,j,k)%background = (13._DP/1000._DP) * thermal(i,j,k)%depth !13K/km for now; must match init_thermal.f90 (this should probably be a global variable)
                end if             
            end do
        end do
    end do

    domain%hmax = hmax
    deallocate(oldtemps)

    return

end subroutine thermal_depth_calculation


