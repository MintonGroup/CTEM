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
e = (rpj/r)**2 - 1.0
delta0 = 12.0*(rpj/r)**2
delta1 = 108.0*(rpj/r)**2
delta  = sqrt(delta1**2 - 4.0*delta0**3)
T = ((delta1 + delta)/2.0)**(1.0/3.0)
S = 0.5*sqrt( 1.0 + (T + delta0/T)/3.0 )
x12 = sqrt(-4.0*S*S + 3.0 + 1.0/S)
x(1) = 0.5 - S + 0.5*x12
x(2) = 0.5 - S - 0.5*x12
x34 = sqrt(-4.0*S*S + 3.0 - 1.0/S)
x(3) = 0.5 + S + 0.5*x34
x(4) = 0.5 + S - 0.5*x34
do i = 1,4
   if (x(i) <= 1.0 .and. x(i) >= 0.0) then
      costheta = x(i)
   end if
end do

!if (present(theta)) then
   z = costheta
!else
!   z = r * (1.0 - costheta) * costheta
!end if

end function regolith_quartic_func
