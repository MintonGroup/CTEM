!**********************************************************************************************************************************
!
!  Unit Name   : thermal_remove
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Removes temperatures of material within the transient cavity that becomes ejecta and moves the temperature distribution 
!                to the surface as in Abramov et al. (2013)
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
subroutine thermal_remove(user,thermal,crater)
    use module_globals
    use module_util
    use module_thermal, EXCEPT_THIS_ONE => thermal_remove
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater

    ! Internal variables
    integer(I4B) :: i, j, k, xpi, ypi, rpix, inc, tdepthpix, reference, maxtdepthpix, npix
    real(DP) :: xp, yp, r, tdepth, maxtdepth, vert, horiz, maxdisp
    type(thermaltype),dimension(:,:,:),allocatable :: oldtemps
    !real(DP),dimension(:,:,:),allocatable :: oldtemps

    allocate(oldtemps,source=thermal(:,:,:))


    rpix = crater%rad / user%pix
    inc = min(rpix,user%gridsize-1)

    maxdisp = 0.06*(crater%fcrat/1000)**1.1 ! eq. 8 in Abramov et al. (2013)

    maxtdepth =  -0.5 * (0.0_DP/crater%rad)**2 + (crater%rad/2._DP)
    maxtdepthpix = maxtdepth / user%zpix

    do j=-inc,inc
        do i=-inc,inc
            ! find distance from crater center
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
   
            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
   
            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! calculate the depth of the transient crater at this pixel, assuming parabolic shape:
            r = sqrt(xp**2 + yp**2)
            tdepth = -0.5 * (r/crater%rad)**2 + (crater%rad/2._DP)
            tdepthpix = tdepth / user%zpix

            horiz = ((r - (0.22*(crater%frad/1000))**2)) / maxdisp !max depth at this horizontal distance

            do k=1,user%zgridsize
                if (thermal(xpi,ypi,k)%depth > 0) then
                    if(thermal(xpi,ypi,k)%depth < maxtdepth) then
                        if (thermal(xpi,ypi,k)%depth < user%zpix) then
                            reference = k
                        end if

                        vert = horiz - ((horiz / (1.25*maxtdepth))*thermal(xpi,ypi,k)%depth) !vertical uplift
                        npix = vert / user%zpix !number of pixels to shift
                        

                        if (k+tdepthpix .gt. user%zgridsize) then !temperature is equal to the background of the deepst voxel
                            thermal(xpi,ypi,k)%temperature = thermal(xpi,ypi,user%zgridsize)%background
                        else
                            thermal(xpi,ypi,k)%temperature = oldtemps(xpi,ypi,k+tdepthpix)%temperature
                        end if
                        
                        ! thermal(xpi,ypi,k)%depth = oldtemps(xpi,ypi,k)%depth
                        ! thermal(xpi,ypi,k)%relative_depth = oldtemps(xpi,ypi,k)%relative_depth
                        ! thermal(xpi,ypi,k)%elevation = oldtemps(xpi,ypi,k)%elevation
                        ! thermal(xpi,ypi,k)%background = oldtemps(xpi,ypi,k)%background
                    end if
                end if
            end do
        end do
    end do

    !call thermal_uplift(user,thermal,crater,oldtemps)

    deallocate(oldtemps)

    return

end subroutine thermal_remove
