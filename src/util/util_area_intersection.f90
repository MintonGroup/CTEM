!**********************************************************************************************************************************
!
!  Unit Name   : util_area_intersection
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the area of intersection between a circle an a square
!  
!
!  Input
!    Arguments : R    = Radius of the circle
!                xbar = x-axis distance of square center from circle center 
!                ybar = y-axis distance of square center from circle center 
!                P    = length of the side of the square
!
!  Output
!    Arguments : area = fraction of the area of the square that intersects circle
!           
! 
!  Notes       : Employs the agorithm of A.D. Groves, Ballistic Research Lab. Memo. Report #1478, April 1963
!
!**********************************************************************************************************************************
function util_area_intersection(R,xbar,ybar,P) result(area)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_area_intersection
implicit none
interface Groves_F
   function Groves_F(U,V) result(F)
   use module_globals
   implicit none
   real(DP),intent(in) :: U,V
   real(DP) :: F
   end function Groves_F
end interface

real(DP),intent(in) :: R,xbar,ybar,P
real(DP) :: area
real(DP) :: a,b,c,d,R2
integer(I4B) :: i

R2 = R**2
area = 0._DP
do i = 1,4
   ! Establish vertices of equivalent first quadrant rectangle
   a = max(0._DP,(-1._DP)**(0.5_DP * (i**2 - i)) * xbar - 0.5_DP * P)
   b = max(0._DP,(-1._DP)**(0.5_DP * (i**2 + i - 2)) * ybar - 0.5_DP * P)
   c = max(0._DP,(-1._DP)**(0.5_DP * (i**2 - i)) * xbar + 0.5_DP * P - a)
   d = max(0._DP,(-1._DP)**(0.5_DP * (i**2 + i - 2)) * ybar + 0.5_DP * P - b)

   ! Evaluate which of the six case we are in
   if ((a**2 + b**2) >= R2) then ! Case I
      cycle
   else if ((a**2 + (b + d)**2 >= R2).and.((a + c)**2 + b**2 >= R2)) then ! Case II
      area = area + 0.5_DP * R2 * Groves_F( a / R, b / R)
   else if (((a + c)**2 + b**2 < R2).and.(a**2 + (b + d)**2 >= R2)) then ! Case III
      area = area + 0.5_DP * R2 * (Groves_F(a / R, b / R) - Groves_F((a + c) / R, b / R))
   else if ((a**2 + (b + d)**2 < R2).and.(((a + c)**2 + b**2) >= R2)) then ! Case IV
      area = area + 0.5_DP * R2 * (Groves_F(a / R, b / R) - Groves_F(a / R, (b + D) / R))
   else if ((a**2 + (b + d)**2 < R2).and.(((a + c)**2 + b**2) < R2).and.((a + c)**2 + (b + d)**2) > R2) then ! Case V
      area = area + 0.5_DP * R2 * (Groves_F(a / R, b / R) - Groves_F((a + c) / R, b / R) - Groves_F(a / R, (b + d) / R))
   else if (((a + c)**2 + (b + d)**2) <= R2) then ! Case VI
      area = area + c * d
   else
      write(*,*) 
      write(*,*) 'Major error in util_area_intersection. This should not happen!'
   end if
end do

area = area / P**2

return
end function util_area_intersection

function Groves_F(U,V) result(F)
use module_globals
implicit none
real(DP),intent(in) :: U,V
real(DP) :: F
real(DP) :: sqU,sqV,UV

sqU = sqrt(1._DP - U**2)
sqV = sqrt(1._DP - V**2)
UV = U*V

F = asin(sqU * sqV - UV) - U * sqU - V * sqV + 2 * UV

return
end function Groves_F

