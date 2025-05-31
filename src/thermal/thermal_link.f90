!**********************************************************************************************************************************
!
!  Unit Name   : thermal_link
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Links temperature from thermal field to thermal history in surflayer
!  
!
!  Input
!    Arguments : 
!
!
!  Output
!    Arguments :
!           
! 
!  Notes       : 
!
!**********************************************************************************************************************************
subroutine thermal_link(user,thermal,surfi)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_link
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(in) :: thermal
    type(surftype),dimension(:,:),intent(inout) :: surfi

    ! Internal variables
    integer(I4B) :: i,j,k,m,regosize,histsize
    real(DP) :: cum_thickness, current_depth, average_depth, current_therm_depth, tplusone, tminusone, temp
    real(DP) :: dplusone, dminusone, dpercent

    ! Executable code
    do i=1,user%gridsize
        do j=1,user%gridsize
            regosize = size(surfi(i,j)%regolayer)
            cum_thickness = 0
            do k=1,regosize
                current_depth = cum_thickness + surfi(i,j)%regolayer(k)%thickness 
                !interpolate between minimum (cum_thickness) and maximum (current_depth) depth for this layer
                average_depth = (cum_thickness + current_depth) / 2.0_DP
                histsize = size(surfi(i,j)%regolayer(k)%thermalhist)
                do m=1,user%zgridsize !Find temperature at the thermal location corresponding to this depth
                    current_therm_depth = thermal(i,j,m)%depth
                    if (current_therm_depth >= cum_thickness) then
                        !Interpolate temperature from thermal gradient
                        if (m > 1) then
                            if (m < user%zgridsize) then
                                temp = thermal(i,j,m)%temperature
                                tplusone = thermal(i,j,m+1)%temperature
                                tminusone = thermal(i,j,m-1)%temperature
                                dplusone = thermal(i,j,m+1)%depth
                                dminusone = thermal(i,j,m-1)%depth
                                if (thermal(i,j,m)%depth > average_depth ) then ! use dplusone, so it's this voxel and dplusone
                                    dpercent = (average_depth - dminusone) / user%zpix
                                    surfi(i,j)%regolayer(k)%thermalhist(histsize)%temperature = dpercent * (temp - tminusone) + temp
                                    exit
                                else
                                    continue
                                end if
                            else
                                surfi(i,j)%regolayer(k)%thermalhist(histsize)%temperature = thermal(i,j,m)%temperature
                            end if
                        else
                            surfi(i,j)%regolayer(k)%thermalhist(histsize)%temperature = thermal(i,j,m)%temperature
                        end if
                    else if (m == user%zgridsize) then !Contingency
                        write(*,*) "Depth of regolayer is greater than full depth of thermal"
                        !Could interpolate based on background temperature and geothermal gradient..?
                        surfi(i,j)%regolayer(k)%thermalhist(histsize)%temperature = ((13._DP/1000._DP) * average_depth ) + thermal(i,j,m)%temperature
                    else
                        continue
                    end if
                end do

                !And age too?
            end do
        end do
    end do
end subroutine thermal_link



