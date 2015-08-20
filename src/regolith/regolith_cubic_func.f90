!**********************************************************************************************************************************
!
!  Unit Name   : regolith_cubic_func
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Cubic function for solving the head's radius of our tangential strema tube
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
function regolith_cubic_func(c1,c2) result(deltar)
use module_globals
use module_regolith, EXCEPT_THIS_ONE => regolith_cubic_func
implicit none
real(DP),intent(in)   :: c1,c2
real(DP)              :: deltar
real(DP)              :: q,r,delta1
real(DP)              :: k1,k2,theta,a,b

q = c1**2/9.0
r = c1**3/27.0 - c2/2.0
delta1 = q**3 - r**2

if (delta1>0._DP) then
   k1 = 27.0_DP*c2 - 2.0_DP*c1**3
   k2 = 2.0_DP*c1**3
   theta = acos(k1/k2)
   deltar = 2.0_DP/3.0_DP*c1*cos(1.0_DP/3.0_DP*theta) - c1/3.0_DP
   !t2 = 2.0_DP/3.0_DP*c1*cos(1.0_DP/3.0_DP*theta - 2.0_DP*PI/3.0_DP) - c1/3.0_DP
   !t3 = 2.0_DP/3.0_DP*c1*cos(1.0_DP/3.0_DP*theta - 4.0_DP*PI/3.0_DP) - c1/3.0_DP
else if (delta1<=0._DP) then
        a = -1.0_DP * sign(1.0_DP,r) * (abs(r) + sqrt(abs(delta1)))**(1.0/3.0)
        b = q/a
        deltar = (a+b) - c1/3.0
        !if (deltar <= 0.0_DP) write(*,*) 'dr=zero',delta1,deltar,a,b,c1,c2
        !if (a == 0.0_DP) write(*,*) 'a=zero',delta1,deltar,a,b,c1,c2
end if

end function regolith_cubic_func
