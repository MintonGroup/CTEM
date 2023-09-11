!**********************************************************************************************************************************
!
!  Unit Name   : util_tscale
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : converts age from true time in Ga to "interval time" based on the lunar Neukum production function (Neukum et al., 2001)
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  This is a Fortran port of the function "T_scale" from craterproduction.py
!
!**********************************************************************************************************************************
function util_tscale(t) result(tscale)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_tscale
    real(DP), intent(in) :: t
    real(DP) :: tscale

    real(DP) :: N1
    real(DP) :: CSFD = 0.0008173348179780709 !this is the result of the shape function at T=1 Ga for the lunar NPF

    N1 = util_npf_timefunc(t)

    tscale = N1 / CSFD
end function util_tscale



