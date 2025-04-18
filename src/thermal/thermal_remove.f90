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
    integer(I4B) :: i, j, k, xpi, ypi, rpix, inc, tdepthpix, reference, maxtdepthpix, npix, k1,k2
    real(DP) :: xp, yp, r, tdepth, maxtdepth, dz, depth_above, z1,z2, T1, T2, frac, old_depth, depth_below
    type(thermaltype),dimension(:,:,:),allocatable :: oldtemps
    !real(DP),dimension(:,:,:),allocatable :: oldtemps

    allocate(oldtemps,source=thermal(:,:,:))


    rpix = crater%rad / user%pix
    inc = min(rpix,user%gridsize-1)

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
                 do k=1,tdepthpix
                    if (thermal(xpi,ypi,k)%depth > 0) then
                        if(thermal(xpi,ypi,k)%depth < maxtdepth) then
                            if (thermal(xpi,ypi,k)%depth < user%zpix) then
                                reference = k
                            end if

                            if (k+tdepthpix .gt. user%zgridsize) then !temperature is equal to the background of the deepst voxel
                                thermal(xpi,ypi,k)%temperature = thermal(xpi,ypi,user%zgridsize)%background
                            else
                                thermal(xpi,ypi,k)%temperature = oldtemps(xpi,ypi,k+tdepthpix)%temperature + oldtemps(xpi,ypi,k+tdepthpix)%warpedbg
                            end if
                            
                            ! thermal(xpi,ypi,k)%depth = oldtemps(xpi,ypi,k)%depth
                            ! thermal(xpi,ypi,k)%relative_depth = oldtemps(xpi,ypi,k)%relative_depth
                            ! thermal(xpi,ypi,k)%elevation = oldtemps(xpi,ypi,k)%elevation
                            ! thermal(xpi,ypi,k)%background = oldtemps(xpi,ypi,k)%background
                        end if
                    end if
                end do
                !do k = 1, user%zgridsize
                !     old_depth = thermal(xpi, ypi, k)%depth
                
                !     ! If voxel is inside the transient crater (material removed), skip it
                !     if (old_depth <= dz) then
                !         thermal(xpi, ypi, k)%temperature = thermal(xpi, ypi, user%zgridsize)%background
                !         cycle
                !     end if
                
                !     ! This voxel is filled by material that used to be deeper — so add dz
                !     depth_below = old_depth + dz
                
                !     ! Search for bracket depths in oldtemps
                !     do k2 = 2, user%zgridsize
                !         if (oldtemps(xpi, ypi, k2)%depth >= depth_below) then
                !             k1 = k2 - 1
                !             exit
                !         end if
                !     end do
                
                !     ! If outside the model, assign background
                !     if (depth_below > oldtemps(xpi, ypi, user%zgridsize)%depth) then
                !         thermal(xpi, ypi, k)%temperature = thermal(xpi, ypi, user%zgridsize)%background
                !     else
                !         z1 = oldtemps(xpi, ypi, k1)%depth
                !         z2 = oldtemps(xpi, ypi, k2)%depth
                
                !         T1 = oldtemps(xpi, ypi, k1)%temperature + oldtemps(xpi, ypi, k1)%warpedbg
                !         T2 = oldtemps(xpi, ypi, k2)%temperature + oldtemps(xpi, ypi, k2)%warpedbg
                
                !         frac = (depth_below - z1) / (z2 - z1)
                
                !         thermal(xpi, ypi, k)%temperature = (1.0_DP - frac) * T1 + frac * T2
                !     end if
                ! end do
            end if
        end do
    end do

    deallocate(oldtemps)

    return

end subroutine thermal_remove
