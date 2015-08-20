!**********************************************************************************************************************************
!
!  Unit Name   : regolith_melt_func
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates volume of melt or vapr 
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
function regolith_melt_func(r,depthb,erad) result(vol)
use module_globals
use module_regolith, EXCEPT_THIS_ONE => regolith_melt_func
implicit none
real(DP),intent(in) :: r,depthb,erad
real(DP) :: vol0, vol
real(DP) :: xints, yints, ymax, theta
real(DP) :: k1,k2,k3 ! for calculating intersected angle

! Shifted melt zone model
vol = 0.0
! cos^2 + k2 * cos + k3 = 0
k1 = 1.0/( 1.0 + 2.0 *depthb/erad )
k2 = -1.0 - k1 
k3 = ( 1.0 + (depthb**2 - r**2)/erad**2 )*k1
theta = acos(-0.5 * k2 - 0.5 * sqrt(k2**2 - 4.0 * k3))
xints = erad * (1.0 - cos(theta)) * cos(theta)
yints = erad * (1.0 - cos(theta)) * sin(theta)
ymax  = erad/4.0 * sqrt(3.0)
vol0 = PI * ( ( r**2 - 0.5 * erad**2) * xints - depthb**2 * xints + depthb * xints**2 + &
       0.5 * erad * xints**2 + 1.0/12.0 * erad**3 *( 1.0 - (1.0 - 4.0*xints/erad)**1.5 ) )

if (yints <= ymax) then
   vol = vol0
else if (yints > ymax) then 
        vol = vol0 + PI/6.0 * erad**3 * (1.0 - 4.0*xints/erad)**1.5  
end if
 
end function regolith_melt_func
