!**********************************************************************************************************************************
!
!  Unit Name   : util_t_from_scale
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : converts age from "interval time" to true time in Ga based on the lunar Neukum production function (Neukum et al., 2001)
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  This is a Fortran port of the function "T_from_scale" from craterproduction.py. It uses the bisect method so is simple but inefficient.
!
!**********************************************************************************************************************************
function util_t_from_scale(scale,start,finish) result(time)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_t_from_scale
    real(DP), intent(in) :: scale, start, finish
    real(DP) :: time

    real(DP) :: tol = 1e-11_DP
    integer(I4B) :: maxiter = 1000
    integer(I4B) :: i
    real(DP) :: a, b, c


    a = start
    b = finish
    time = -1.
    i = 0
    temp = 0._DP

    if (abs(temp-scale)<tol) then
        time = 0.0_DP
    else
        do while(i .lt. maxiter)
            c = (a+b)/2
            temp = util_tscale(c)

            if (abs(temp-scale)<tol) then
                time = c
                exit
            else if ((temp-scale)>0) then
                b = c
                i = i + 1
            else if ((temp-scale)<0) then
                a = c
                i = i + 1
            else if (abs(a-b)<tol) then
                write(*,*) "ERROR in util_t_from_scale: Convergence failed!"
                write(*,*) scale, start, finish
                exit
            end if
        end do
    end if

    if (time .lt. 0) then
        write(*,*) "ERROR in util_t_from_scale: Maximum iterations reached!"
        write(*,*) scale, maxiter
    end if
end function util_t_from_scale


