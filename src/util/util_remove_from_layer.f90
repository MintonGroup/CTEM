!**********************************************************************************************************************************
!
!  Unit Name   : util_remove_from_layer
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Removes crater from layer
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
subroutine util_remove_from_layer(surfi,layer)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_remove_from_layer
implicit none

! Arguments
type(surftype),intent(inout) :: surfi
integer(I4B),intent(in) :: layer

! Internals

surfi%diam(layer) = 0.0_DP
surfi%xl(layer) = 0.0_SP
surfi%yl(layer) = 0.0_SP

end subroutine util_remove_from_layer
