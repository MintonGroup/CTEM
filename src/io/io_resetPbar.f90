!**********************************************************************************************************************************
!
!  Unit Name   : io_resetPbar
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : places information on the progress bar
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       : 
!
!**********************************************************************************************************************************
subroutine io_resetPbar()
use module_globals
use module_io, EXCEPT_THIS_ONE => io_resetPbar 
implicit none
integer(I4B) :: k

pbarpos=0
do k=1,pbarsize
   pbarchar(k:k)="="
end do

return 
end subroutine io_resetPbar

