!**********************************************************************************************************************************
!
!  Unit Name   : module_porosity
!  Unit Type   : module
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for porosity creation
!
!  Notes       :  
!
!  Developer   : Toshi Hirabayashi
!**********************************************************************************************************************************
module module_porosity
use module_globals
implicit none
public 
save

   interface
      subroutine porosity_form_interior(user, surfi, crater, lradsq)
      use module_globals
      implicit none
      type(usertype),intent(in)    :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in)  :: crater
      real(DP),intent(in)          :: lradsq
      end subroutine porosity_form_interior
   end interface

end module
