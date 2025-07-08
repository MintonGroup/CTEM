!**********************************************************************************************************************************
!
!  Unit Name   : thermal_ejecta_average.f90
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Averages non-vapor temperatures within the transient crater region of the thermal structure
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : 
!           
! 
!  Notes       :  This is what Abramov et al. (2013) do, but might not be the best approach here. This is a first-pass exercise
!
!**********************************************************************************************************************************
subroutine thermal_ejecta_average(user,surf,crater,thermal,avgtemp)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_ejecta_average
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(surftype),dimension(:,:),intent(in) :: surf
    type(cratertype),intent(in) :: crater
    type(thermaltype),dimension(:,:,:),intent(in) :: thermal
    real(DP),intent(out) :: avgtemp

    ! Internal variables
    integer(I4B) :: i,j,k,inc,rpix,xpi,ypi,nvox,m,n,nz
    real(DP) :: xp,yp,maxtdepth,excdepth,tdepth,r,depth
    real(DP),dimension(:),allocatable :: temps,nonzero

    ! Executable code
    n = 0

    maxtdepth = 0.25* (2*crater%rad)
    excdepth = 0.5 * maxtdepth

    !Find the radial extent of the transient crater in pixels
    rpix = crater%rad / user%pix
    inc = min(rpix,user%gridsize-1)

    ! First pass to find number of voxels to average

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

            if (tdepth > excdepth) then
                depth = excdepth
            else
                depth = tdepth
            end if

            nvox = int(depth/user%zpix)
            n = n + nvox
        end do
    end do

    allocate(temps(n))
    m = 1

    temps(:) = 0.0_DP

    ! Second pass to actually get the temperature value of each voxel that is not vapor

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

            if (tdepth > excdepth) then
                depth = excdepth
            else
                depth = tdepth
            end if

            nvox = int(depth/user%zpix)

            do k=1,nvox
                if (thermal(xpi,ypi,k)%temperature < 3327.0) then
                    temps(m) = thermal(xpi,ypi,k)%temperature
                else
                    temps(m) = 0.0_DP
                end if
                m = m + 1
            end do
        end do
    end do

    ! find the number of non-zero values in the array

    nz = 0

    do i=1,n
        if (temps(i) > 0.0_DP) nz = nz + 1
    end do

    allocate(nonzero(nz))

    m = 1

    do i=1,n
        if (temps(i) > 0.0_DP) then
            nonzero(m) = temps(i)
            m = m + 1
        end if
    end do

    ! Finally get the average
    if (size(nonzero)>1) then
        avgtemp = sum(nonzero) / size(nonzero)
    else
        avgtemp = 0.0
    end if

    deallocate(temps,nonzero)

    return
end subroutine thermal_ejecta_average




    