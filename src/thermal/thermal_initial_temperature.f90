!**********************************************************************************************************************************
!
!  Unit Name   : thermal_initial_temperature
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Add the initial thermal distribution from an impact based on the scaling law found in Abramov et al. (2013)
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
subroutine thermal_initial_temperature(user,crater,thermi,distance)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_initial_temperature
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(cratertype),intent(in) :: crater
    type(thermaltype),intent(inout) :: thermi
    real(DP),intent(in) :: distance

    ! Internal variables
    real(DP) :: term1,term2,term3, deltaT, deltaEw, P, k, A, V0, K0, n, C, Cprime, oldtemp

    ! Executable code
    V0 = 1.0_DP / user%trho_r
    K0 = 19.3e9 !47.5e9 !zero-pressure bulk modulus for anorthosite (data from McQueen, 1967, method from Kieffer and Simmonds, 1980) <--from HEATING (Abramov et al., 2013)
    n = 5.5 !3.49 !constant for anorthosite (data from McQueen, 1967, method from Kieffer and Simmonds, 1980)  <--from HEATING (Abramov et al., 2013)

    k = 0.625_DP*log10(crater%impvel/1000._DP) + 1.25 !Equation 3 in Abramov et al. (2013)
    A = 0.25_DP * user%prho * crater%impvel**2 * crater%sinimpang ! Equation 4 in Abramov et al. (2013)
                                                                  !prho assumed to be the same as target density
                                                                  ! Collins et al. (2002) may have the derivation for this equation
                                                                  ! (in case I need to modify it for when densities are different)
    C = 800. !820 is value for anorthosite in HEATING (Jones, 2015). 800 for bssalt
    P = A*(distance/crater%imprad)**(-k) !Equation 2 in Abramov et al. (2013)

    term1 = 0.5_DP * (P * V0 - (2 * K0 * V0) / n)
    term2 = 1 - ((P * n / K0) + 1)**(-1/n)
    term3 = (K0 * V0 / (n * (1 - n))) * (1 - ((P * n / K0) + 1)**(1 - (1/n)))
    deltaEw = term1 * term2 + term3 !Equation 1 in Abramov et al. (2013)
    deltaT = deltaEw / C 
    oldtemp = thermi%temperature
    thermi%temperature = thermi%temperature + deltaT

    if (thermi%temperature >= user%tsolidus .and. thermi%temperature <= user%tliquidus) then !Add latent heat calculation
        Cprime = C + (330. / (user%tliquidus - user%tsolidus)) ! 330 is latent heat of fusion for basalt as used by Melosh 2000
        deltaT = deltaEw / Cprime
        thermi%temperature = oldtemp + deltaT
    end if

    return
end subroutine thermal_initial_temperature