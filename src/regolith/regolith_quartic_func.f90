!**********************************************************************************************************************************
!
!  Unit Name   : regolith_quartic_func
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Quartic formula for depth in a stream tube at a specific subpixel
!  
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
function regolith_quartic_func(rpj,r) result(z)
use module_globals
use module_regolith, EXCEPT_THIS_ONE => regolith_quartic_func
implicit none
real(DP),intent(in)              :: rpj,r
!real(DP),intent(inout),optional  :: theta
real(DP)                         :: z
real(DP),parameter    :: a = 1.0
real(DP),parameter    :: b = -2.0
real(DP),parameter    :: c = 0.0
real(DP),parameter    :: d = 2.0
real(DP),parameter    :: p = -1.5
real(DP),parameter    :: q = 1.0
real(DP)              :: e,delta0,delta1,delta,S,T,x12,x34,costheta
real(DP),dimension(4) :: x
integer(I4B)          :: i

costheta = 0.0
e = (rpj/r)**2 - 1.0_DP
delta0 = 12.0_DP*(rpj/r)**2
delta1 = 108.0_DP*(rpj/r)**2
delta  = sqrt(delta1**2 - 4.0*delta0**3)
T = ((delta1 + delta)/2.0_DP)**(1.0_DP/3.0_DP)
S = 0.5*sqrt( 1.0_DP + (T + delta0/T)/3.0_DP )
x12 = sqrt(-4.0_DP*S*S + 3.0_DP + 1.0_DP/S)
x(1) = 0.5 - S + 0.5*x12
x(2) = 0.5 - S - 0.5*x12
x34 = max(-4*S**2 + 3.0_DP - 1.0_DP/S,0.0_DP)
x34 = sqrt(x34)
x(3) = 0.5_DP + S + 0.5_DP*x34
x(4) = 0.5_DP + S - 0.5_DP*x34
do i = 1,4
   if (x(i) <= 1.0_DP .and. x(i) >= 0.0_DP) then
      costheta = x(i)
   end if
end do

!if (present(theta)) then
z = costheta
!else
!   z = r * (1.0 - costheta) * costheta
!end if

end function regolith_quartic_func
