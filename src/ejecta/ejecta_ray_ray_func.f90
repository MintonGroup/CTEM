!****f* ejecta/ejecta_ray_ray_func
! Name
!   ejecta_ray_ray_func -- Calculate ejecta ray
! SYNOPSIS
!   This uses
!   * module_globals
!   * module_ejecta
!
!   ans = ejecta_ray_ray_func()
!
! DESCRIPTION
!   
!   This function was separated into an independent function for debugging purposes.
!  
! ARGUMENTS
!   Input
!   * theta      -- 
!   * thetar     -- 
!   * r    -- 
!   * n        -- 
!   * w      -- 
!
!   Output
!   * ans  -- THe result of this function, used in ejecta_ray_pattern.f90
! 
!***

!**********************************************************************************************************
function ejecta_ray_ray_func(theta,thetar,r,n,w) result(ans)
    use module_globals
    use module_ejecta, EXCEPT_THIS_ONE => ejecta_ray_ray_func
    implicit none
    real(DP) :: ans
    real(DP),intent(in) :: theta,thetar,r,w
    integer(I4B),intent(in) :: n
    real(DP) :: thetap,thetapp,a,b,c,dtheta

    c = w / r
    b = thetar 
    dtheta = min(2*pi - abs(theta - b),abs(theta - b))
    a = sqrt(2 * pi) / (n * c * erf(pi / (2 *sqrt(2._DP) * c)))
    ans = a * exp(-dtheta**2 / (2 * c**2))

!return
end function ejecta_ray_ray_func