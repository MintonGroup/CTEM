!**********************************************************************************************************************************
!
!  Unit Name   : util_poisson
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Draws a random number from a Poisson distribution
!  
!
!  Input
!    Arguments : mu  - The mean of the distribution
!                first  - Set to true if this is the first time a deviate is calculated from a given mu
!
!  Output
!    Arguments : ival - The Poisson random deviate
!           
! 
!  Notes       :
!
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!     Adapted for CTEM by David Minton
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate POIsson random deviate

!                            Function

! Generates a single random deviate from a Poisson distribution with mean lambda.

!**********************************************************************************************************************************
function util_poisson(mu,poisson_first) result(ival)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_poisson
implicit none
interface random_normal
   function random_normal() result(fn_val)
   use module_globals
   implicit none
   real(DP) ::fn_val
   end function random_normal
end interface
interface random_exponential
   function random_exponential() result(fn_val)
   use module_globals
   implicit none
   real(DP)  :: fn_val
   end function random_exponential
end interface


!     .. Scalar Arguments ..
real(DP), intent(in)    :: mu
logical, intent(in),optional :: poisson_first
integer(I8B)             :: ival
!     ..
!     .. Local Scalars ..
real(DP)          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                 omega, px, py, t, u, v, x, xx
real(DP), save    :: s, d, p, q, p0
integer(I8B)       :: j, k, kflag
logical, save :: full_init
integer(I8B), save :: l, m
!     ..
!     .. Local Arrays ..
real(DP), save    :: pp(35)
!     ..
!     .. Data statements ..
real(DP), parameter :: a0 = -.5_DP, a1 = .3333333_DP, a2 = -.2500068_DP, a3 = .2000118_DP,  &
                   a4 = -.1661269_DP, a5 = .1421878_DP, a6 = -.1384794_DP,   &
                   a7 = .1250060_DP

real(DP), parameter :: fact(10) = (/ 1._DP, 1._DP, 2._DP, 6._DP, 24._DP, 120._DP, 720._DP, 5040._DP,  &
                                 40320._DP, 362880._DP /)
real(DP), parameter  :: zero = 0.0_DP, half = 0.5_DP, one = 1.0_DP, two = 2.0_DP
logical :: first
                      

!     ..
!     .. Executable Statements ..
if (present(poisson_first)) then
   first = poisson_first
else
   first = .true.
end if

if (mu > 10.0_DP) then
!     C A S E  A. (RECALCULATION OF S, D, L if MU HAS CHANGED)

  if (first) then
    s = sqrt(mu)
    d = 6.0_DP * mu**2

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=ifIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

    l = mu - 1.1484_DP
    full_init = .false.
  end if


!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

  g = mu + s * random_normal()
  if (g > 0.0_DP) then
    ival = g

!     STEP I. IMMEDIATE ACCEPTANCE if ival IS LARGE ENOUGH

    if (ival>=l) return

!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

    fk = ival
    difmuk = mu - fk
    call random_number(u)
    if (d * u >= difmuk**3) return
  end if

!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF parameterS if NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

  if (.not. full_init) then
    omega = .3989423_DP / s
    b1 = .4166667e-1_DP / mu
    b2 = .3_DP * b1 * b1
    c3 = .1428571_DP * b1 * b2
    c2 = b2 - 15._DP  *c3
    c1 = b1 - 6._DP * b2 + 45._DP * c3
    c0 = 1._DP - b1 + 3._DP * b2 - 15._DP * c3
    c = .1069_DP / mu
    full_init = .true.
  end if

  if (g < 0.0_DP) goto 50

!             'SUBROUTINE' F IS callED (KFLAG=0 FOR CORRECT return)

  kflag = 0
  goto 70

!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

  40 if (fy-u*fy <= py*exp(px-fx)) return

!     STEP E. expONENTIAL SAMPLE - random_exponential() FOR STANDARD expONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (if T <= -.6744 then PK < FK FOR ALL MU >= 10.)

  50 e = random_exponential()
  call random_number(u)
  u = u + u - one
  t = 1.8_DP + sign(e, u)
  if (t <= (-.6744_DP)) goto 50
  ival = mu + s * t
  fk = ival
  difmuk = mu - fk

!             'SUBROUTINE' F IS callED (KFLAG=1 FOR CORRECT return)

  kflag = 1
  goto 70

!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

  60 if (c*abs(u) > py*exp(px+e) - fy*exp(fx+e)) goto 50
  return

!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
!             CASE ival < 10 USES FACTORIALS FROM TABLE FACT

  70 if (ival>=10) goto 80
  px = -mu
  py = mu**ival / fact(ival+1)
  goto 110

!             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
!             A0-A7 FOR ACCURACY WHEN ADVISABLE
!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

  80 del = .8333333E-1_DP / fk
  del = del - 4.8_DP * del**3
  v = difmuk/fk
  if (abs(v)>0.25_DP) then
    px = fk*log(one + v) - difmuk - del
  else
    px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  end if
  py = .3989423_DP / sqrt(fk)
  110 x = (half - difmuk)/s
  xx = x*x
  fx = -half*xx
  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
  if (kflag <= 0) goto 40
  goto 60

!---------------------------------------------------------------------------
!     C A S E  B.    mu < 10
!     START NEW TABLE AND CALCULATE P0 if NECESSARY

else
  if (first) then
    m = max(1, int(mu))
    l = 0
    p = exp(-mu)
    q = p
    p0 = p
  end if

!     STEP U. UNifORM SAMPLE FOR INVERSION METHOD

  do
    call random_number(u)
    ival = 0
    if (u <= p0) return

!     STEP T. TABLE COMPARISON UNTIL THE end PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458=PP(9) FOR MU=10)

    if (l == 0) goto 150
    j = 1
    if (u > 0.458_DP) j = MIN(l, m)
    do k = j, l
      if (u <= pp(k)) goto 180
    end do
    if (l == 35) cycle

!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q=PP(K)

    150 l = l + 1
    do k = l, 35
      p = p*mu / k
      q = q + p
      pp(k) = q
      if (u <= q) goto 170
    end do
    l = 35
  end do

  170 l = k
  180 ival = k
  return
end if

return
end function util_poisson

function random_normal() result(fn_val)
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.
use module_globals
implicit none
real(DP) :: fn_val
!     Local variables
real(DP) :: s = 0.449871_DP, t = -0.386595_DP, a = 0.19600_DP, b = 0.25472_DP,    &
            r1 = 0.27597_DP, r2 = 0.27846_DP, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

do
   call random_number(u)
   call random_number(v)
   v = 1.7156_DP * (v - 0.5_DP)

   !     Evaluate the quadratic form
   x = u - s
   y = abs(v) - t
   q = x**2 + y * (a * y - b * x)

   !     Accept P if inside inner ellipse
   if (q < r1) exit
   !     Reject P if outside outer ellipse
   if (q > r2) cycle
   !     Reject P if outside acceptance region
   if (v**2 < -4.0_DP * log(u) * u**2) exit
end do

!     Return ratio of P's coordinates as the normal deviate
fn_val = v / u

return
end function random_normal

function random_exponential() result(fn_val)
use module_globals
implicit none
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

real(DP)  :: fn_val

!     Local variable
real(DP)  :: r

do
  call random_number(r)
  if (r > 0.0_DP) exit
end do

fn_val = -log(r)
return
end function random_exponential


