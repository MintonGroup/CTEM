!**********************************************************************************************************************************
!
!  Unit Name   : regolith_circle_sector_func
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : circle sector approximation for the head part of a stream tube
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
function regolith_circle_sector_func(deltar,zstart,zend) result(vhead)
use module_globals
use module_regolith, EXCEPT_THIS_ONE => regolith_circle_sector_func
implicit none
real(DP),intent(in) :: deltar,zstart,zend 
real(DP)            :: vhead
real(DP)            :: zi,zf

zi = abs(1.0_DP - zstart/deltar)
zf = abs(1.0_DP - zend/deltar)

vhead = 1.0_DP/3.0_DP*((1.0_DP-zf**2)**(1.5)-(1.0_DP-zi**2)**(1.5)) &
        + (zf*acos(zf) - zi*acos(zi)+sqrt(1.0_DP-zi**2)-sqrt(1.0_DP-zf**2))
!write(*,*) PI*deltar**3,vhead
if (zstart<=deltar .and. zend<=deltar) then 
   vhead = PI*(zi-zf) + vhead
   !write(*,*) 'abo',zi,zf,zstart,zend,deltar,vhead
else if (zstart>deltar .and. zend>deltar) then     
        vhead = vhead
        !write(*,*) 'bel',zi,zf,zstart,zend,deltar,vhead
else 
    vhead = PI*zi + vhead
    !write(*,*) 'ibe',zi,zf,zstart,zend,deltar,vhead
end if

vhead = deltar**3 * vhead

end function regolith_circle_sector_func
