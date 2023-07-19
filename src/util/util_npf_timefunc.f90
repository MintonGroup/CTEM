!**********************************************************************************************************************************
!
!  Unit Name   : util_npf_timefunc
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : The time function of the lunar Neukum production function (Neukum et al., 2001)
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :
!
!**********************************************************************************************************************************
function util_npf_timefunc(T) result(N1)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_npf_timefunc
    real(DP), intent(in) :: T
    real(DP) :: N1

    N1 = 5.44e-14 * (exp(6.93*T)-1) + 8.17e-4*T
end function util_npf_timefunc
