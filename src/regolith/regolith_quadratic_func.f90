!**********************************************************************************************************************************
!
!  Unit Name   : regolith_Quadratic_func
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
function regolith_quadratic_func(z,erad,ri,rip1,rstart) result(r)
use module_globals
use module_regolith, EXCEPT_THIS_ONE => regolith_quadratic_func
implicit none
real(DP),intent(in)              :: z,erad,ri,rip1,rstart
!real(DP),intent(in),optional     :: theta
real(DP)           :: r

!Internal Arguments
real(DP) :: theta1,theta2,r1,r2

r = rstart

!if (z<=erad/4.0_DP) then

theta1 = acos(0.5_DP + 0.5_DP*sqrt(1.0_DP - 4.0_DP*z/erad))
theta2 = acos(0.5_DP - 0.5_DP*sqrt(1.0_DP - 4.0_DP*z/erad))

r1 = erad * (1.0_DP - cos(theta1)) * sin(theta1)
r2 = erad * (1.0_DP - cos(theta2)) * sin(theta2)

if (r1>=ri .and. r1<=rip1) then
   r = r1
else if (r2>=ri .and. r2<=rip1) then
   r = r2
end if

!end if

end function regolith_quadratic_func
