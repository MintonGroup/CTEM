!**********************************************************************************************************************************
!
!  Unit Name   : util_add_to_layer
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Adds elevation change to layer
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
subroutine util_add_to_layer(user,surfi,fcrat,xl,yl)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_add_to_layer
implicit none

! Arguments
type(usertype),intent(in) :: user
type(surftype),intent(inout) :: surfi
real(DP),intent(in) :: fcrat
real(SP),intent(in) :: xl,yl

! Internals
integer(I4B) :: layer,l

! Executable code
layer=0
!  add to successive layers, destroying smaller craters if they exist
do l=user%numlayers,1,-1
   ! mark layer as available if nothing as big is now recorded there:
   if ((fcrat * COOKIESIZE) > surfi%diam(l) ) then 
      call util_remove_from_layer(surfi,l)
      layer = l
   end if
end do

! Emplace a new crater if requested and a spot is available
if (layer>0) then
   surfi%diam(layer) = fcrat
   surfi%xl(layer) = xl
   surfi%yl(layer) = yl
else
   write(*,*)
   write(*,*) 'WARNING! No free layer to add crater pixel. Consider increasing NUMLAYERS'
   !write(*,*) 'Crater size in pixels: ',int(fcrat/user%pix)
end if

end subroutine util_add_to_layer
