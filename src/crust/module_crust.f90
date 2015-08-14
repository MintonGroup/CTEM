!**********************************************************************************************************************************
!
!  Unit Name   : module_crust
!  Unit Type   : module
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for crust creation
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_crust
use module_globals
implicit none
public 
save

   interface
      subroutine crust_thin(user,surf,crater,domain,mdepth)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: mdepth
      end subroutine crust_thin
   end interface

   interface
      subroutine crust_cta_generate(user,surfi,crater,domain,lradsq,depth,mdepth,trans_depth)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: lradsq
      real(DP),intent(in) :: depth
      real(DP),intent(in) :: mdepth
      real(DP),intent(in) :: trans_depth
      end subroutine crust_cta_generate
   end interface

end module
