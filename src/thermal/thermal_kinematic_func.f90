!**********************************************************************************************************************************
!
!  Unit Name   : thermal_kinematic_func.f90
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds the flight time for a parcel of ejecta and computes Ar loss based on assumed heating
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       : Angle in degrees
!
!**********************************************************************************************************************************
function thermal_kinematic_func(vesq,angle,lrad) result(loss)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_kinematic_func
    implicit none

    real(DP),intent(in) :: vesq,angle,lrad
    real(DP) :: loss

    !Internal
    real(DP) :: t, v0, theta, temp, f, dr2

    !Executable Code

    ! Find time
    v0 = sqrt(vesq)
    t = lrad / (v0 * cosd(angle))

    ! Use this time to calculate Ar loss from a specific temperature

    temp = 1000 !Not quite melt

    dr2  = exp(-2.10_DP*(1e4_DP/(temp+223))+8.05_DP) ! Assumes -50C for surface temperature; change to +273 for 0C
    f = ((6.0_DP/PI**(1.5_DP))*(((PI)**2.0_DP)*dr2*t)**(0.5_DP))-(((3/(PI**2.0_DP))*((PI)**2.0_DP)*dr2*t))
    if (f < 0._DP .or. f > 0.85_DP) then
        f = 1.0_DP-(6.0_DP/PI**2.0_DP)*exp((-(PI)**2.0_DP)*dr2*t)
    end if
    if (f >= 0._DP .and. f <= 1._DP) then
        loss = f
    else
        loss = -1.0_DP
    end if
end function thermal_kinematic_func



