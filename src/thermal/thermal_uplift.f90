!**********************************************************************************************************************************
!
!  Unit Name   : thermal_uplift
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Adds uplift to central region following Abramov et al. (2013) p. 233-234
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
subroutine thermal_uplift(user,thermal,crater,oldtemps)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_uplift
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater
    type(thermaltype),dimension(:,:,:), intent(in) :: oldtemps

    ! Internal variables

    integer(I4B) :: i, j, k, xpi, ypi, rpix, inc, tdepthpix, reference, maxtdepthpix, npix
    real(DP) :: xp, yp, r, tdepth, maxtdepth, maxdisp, vert, horiz

    maxtdepth =  -0.5 * (0.0_DP/crater%rad)**2 + (crater%rad/2._DP)
    maxtdepthpix = maxtdepth / user%zpix

    maxdisp = 0.06*(crater%fcrat/1000)**1.1 ! eq. 8 in Abramov et al. (2013)

    rpix = (0.22*crater%frad) / user%pix !lateral extent of the uplift
    inc = min(rpix,user%gridsize-1)

    
