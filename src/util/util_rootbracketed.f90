!**********************************************************************************************************************************
!
!  Unit Name   : util_rootbracketed
!  Unit Type   : function (logical)
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description :  TRUE if x1*x2 negative
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
function util_rootbracketed(x1,x2) result (resultat)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_rootbracketed
implicit none

! Arguments
real(DP),intent(in) :: x1,x2
logical :: resultat

! Executable code

  if ((x1 > 0.and.x2 > 0).or.(x1 < 0.and.x2 < 0)) then 
    resultat = .false.
  else
    resultat = .true.
  endif

  return
end function util_rootbracketed

