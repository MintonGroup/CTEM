!****f* ejecta/ejecta_ray_pattern_func
! Name
!   ejecta_ray_pattern_func -- Calculate ejecta ray pattern
! SYNOPSIS
!   This uses
!   * module_globals
!   * module_ejecta
!
!   ans = ejecta_ray_pattern_func()
!
! DESCRIPTION
!   
!   This function was separated into an independent function for debugging purposes.
!  
! ARGUMENTS
!   Input
!   * theta      -- 
!   * r     -- 
!   * rmin    -- 
!   * rmax        -- 
!   * thetari      -- 
!   * ej     -- 
!
!   Output
!   * ans  -- THe result of this function, used in ejecta_ray_pattern.f90
! 
!***

!**********************************************************************************************************
function ejecta_ray_pattern_func(theta,r,rmin,rmax,thetari,ej) result(ans)
    use module_globals
    use module_ejecta, EXCEPT_THIS_ONE => ejecta_ray_pattern_func
    implicit none
    real(DP) :: ans
    real(DP),intent(in) :: r,rmin,rmax,theta
    real(DP),dimension(:),intent(in) :: thetari
    logical,intent(in) :: ej
    real(DP) :: a,c
    real(DP) :: thetar,rw,rw0,rw1
    real(DP) :: f,rtrans,length,rpeak,minray,FF
    integer(I4B) :: n,i
 
 
    minray = rmin * 3
 
    if (r > rmax) then
       ans = 0._DP
    else if (r < 1.0_DP) then
       if (ej) then
          ans = 1.0_DP
       else
          ans = 0.0_DP
       end if
    else
       rw0 = rmin * pi / Nraymax / 2
       rw1 = 2 * pi / Nraymax
       rw = rw0 * (1._DP - (1.0_DP - rw1 / rw0) * exp(1._DP - (r / rmin)**2)) ! equation 40 Minton et al. 2019
       n = max(min(floor((Nraymax**rayp - (Nraymax**rayp - 1) * log(r/minray) / log(rray/minray))**(1._DP/rayp)),Nraymax),1) ! Exponential decay of ray number with distance
       ans = 0._DP
       rtrans = r - 1.0_DP
       c = rw / r
       a = sqrt(2 * pi) / (n * c * erf(pi / (2 *sqrt(2._DP) * c))) !equation 39 Minton et al., 2019
       do i = 1,Nraymax
          length = minray * exp(log(rray/minray) * ((Nraymax - i + 1)**rayp - 1_DP) / ((Nraymax**rayp - 1)))
          rpeak = (length - 1_DP) * 0.5_DP
          if (ej) then
             FF = 1.0_DP
             if (r > length) then
                f = 0.0_DP
             else
                f = a 
             end if
          else
             FF = rayfmult * (20 / rmax)**(0.5_DP) * 0.25_DP 
             f = FF * fpeak * (rtrans / rpeak)**rayq * exp(1._DP / rayq * (1.0_DP - (rtrans / rpeak)**rayq)) !equation 42 Minton et al. 2019
          end if
          ans = ans + ejecta_ray_ray_func(theta,thetari(i),r,n,rw) * f / a 
       end do
    end if 
 
end function ejecta_ray_pattern_func