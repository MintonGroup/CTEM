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
subroutine thermal_remove(user,thermal,crater,prev)
    use module_globals
    use module_util
    use module_thermal, EXCEPT_THIS_ONE => thermal_remove
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater
    real(DP),dimension(:,:,:),allocatable,intent(in) :: prev

    ! Internal variables
    integer(I4B) :: i, j, k, xpi, ypi, rpix, inc, tdepthpix, reference, maxtdepthpix, npix, k1,k2
    real(DP) :: xp, yp, r, tdepth, maxtdepth, dz, depth_above, z1,z2, T1, T2, frac, old_depth, depth_below, limit
    type(thermaltype),dimension(:,:,:),allocatable :: oldtemps
    !real(DP),dimension(:,:,:),allocatable :: oldtemps

    allocate(oldtemps,source=thermal(:,:,:))


    rpix = crater%rad / user%pix
    inc = min(rpix,user%gridsize-1)
    limit = 13.0*crater%imprad

    maxtdepth = crater%rad * ((-0.5 * (0.0_DP/crater%rad)**2) + 0.5) ! Parabolic relationship between r and depth of transient crater: 
                                                                    ! depth h = -0.5(r/R)**2 + (1/2) where r is radial distance and R is radius. 
                                                                    ! Max depth of transient crater is assumed to be at r=0, and goes to 0 at r=R.
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
            r = sqrt((crater%xl-xp)**2 + (crater%yl-yp)**2)
            tdepth = -0.5 * (r/crater%rad)**2 + 0.5 ! General use case of parabolic relationship described above

            dz = tdepth * crater%rad ! continuous uplift

            tdepthpix = dz / user%zpix

            if (dz > 0.0_DP) then
                 do k=1,user%zgridsize
                    if (thermal(xpi,ypi,k)%depth > 0) then
                        if(thermal(xpi,ypi,k)%depth < maxtdepth) then
                            if (thermal(xpi,ypi,k)%depth < user%zpix) then
                                reference = k
                            end if
                            if (k+tdepthpix .gt. user%zgridsize) then !temperature is equal to the background of the deepst voxel
                                thermal(xpi,ypi,k)%temperature = thermal(xpi,ypi,user%zgridsize)%background
                            else
                                if (thermal(xpi,ypi,k)%depth < limit) then !Remove transient stuff and shift
                                    thermal(xpi,ypi,k)%temperature = oldtemps(xpi,ypi,k+tdepthpix)%temperature + oldtemps(xpi,ypi,k+tdepthpix)%warpedbg
                                else !if k > limit, make thermal the background (for now, it should actually be the value of "prev")
                                    thermal(xpi,ypi,k)%temperature = prev(xpi,ypi,k)
                                end if
                            end if
                        else
                            thermal(xpi,ypi,k)%temperature = prev(xpi,ypi,k)
                        end if
                    end if
                end do
            end if
        end do
    end do

    deallocate(oldtemps)

    ! Change the warpedbg value back to 0
    ! do i=1,user%gridsize
    !     do j=1,user%gridsize
    !         do k=1,user%zgridsize
    !             thermal(i,j,k)%warpedbg = 0
    !         end do
    !     end do
    ! end do

    thermal(:,:,:)%warpedbg = 0.0_DP

    return

end subroutine thermal_remove
