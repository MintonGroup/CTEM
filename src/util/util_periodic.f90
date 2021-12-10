!**********************************************************************************************************************************
!
!  Unit Name   : util_periodic
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Implements periodic boundary conditions
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
pure subroutine util_periodic(x,y,side)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_periodic
implicit none

! Arguments
integer(I4B),intent(inout) :: x,y
integer(I4B),intent(in) :: side

! Executable code

do while (x < 1)
   x = x + side
end do
do while (x > side )
   x = x - side
end do

do while (y < 1)
   y = y + side
end do
do while (y > side )
   y = y - side
end do

return
end subroutine util_periodic
