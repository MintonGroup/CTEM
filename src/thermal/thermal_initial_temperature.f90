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
subroutine thermal_initial_temperature(user,crater,thermal)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_initial_temperature
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(cratertype),intent(in) :: crater
    type(thermaltype),intent(inout) :: thermal