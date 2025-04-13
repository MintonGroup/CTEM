!**********************************************************************************************************************************
!
!  Unit Name   : thermal_warp
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Warps background temperature (bottom left F.3, Abramov et al. 2013)
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
subroutine thermal_warp(user,thermal,crater)
    use module_globals
    use module_util
    use module_thermal, EXCEPT_THIS_ONE => thermal_warp
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater

    ! Internal variables
    integer(I4B) :: i,j,k,inc,xpi,ypi, maxzpix,rpix!,dummy
    real(DP) :: Rcp, trans_depth, maxdisp, r, uz, maxz, z, xp, yp

    ! Executable code

    trans_depth = 0.5 * crater%rad ! 0.25 * 2 * crater%rad
    Rcp = 0.22 * crater%frad ! lateral extent of uplift, p.8 A13
    maxdisp = (0.06*(crater%fcrat/1000)**1.1) * 1000 ! eq. 8 in Abramov et al. (2013), converted to m

    rpix = crater%rad / user%pix
    inc = min(rpix,user%gridsize-1)


    do j = -inc,inc
        do i =-inc,inc
            ! find distance from crater center
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
   
            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
   
            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            r = sqrt((crater%xl-xp)**2 + (crater%yl-yp)**2)

            maxz = 1.25*trans_depth
            maxzpix = maxz / user%zpix

            if (r <= Rcp) then
                do k=1,maxzpix
                    z = k * user%zpix
                    uz = maxdisp * (1.0_DP - (z/maxz)) * (1.0_DP - (r/Rcp)**2) ! vertical displacement
                    uz = min(z,uz) !No negative values-- things above the surface are removed
                    thermal(xpi,ypi,k)%warp = z - uz
                    !dummy = 1 ! for debugging the line above
                end do
            end if
        end do
    end do

    ! Assign background temperature from warping

    ! Make warp go back to 0

    ! do i=1,user%gridsize
    !     do j=1,user%gridsize
    !         do k=1,user%zgridsize
    !             thermal(i,j,k)%warp = 0.0_DP
    !         end do
    !     end do
    ! end do
    return
end subroutine thermal_warp










































