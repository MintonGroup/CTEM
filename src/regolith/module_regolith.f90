!**********************************************************************************************************************************
!
!  Unit Name   : module_crater
!  Unit Type   : module
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for regolith layer !  tracking
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_regolith
use module_globals
implicit none
public 
save

   interface
      subroutine regolith_transport(user,surf,crater,domain,ejb,ejtble,xp,yp)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      type(ejbtype),dimension(:),intent(in)    :: ejb
      integer(I4B),intent(in) :: ejtble
      real(DP),intent(in) :: xp,yp
      end subroutine regolith_transport
   end interface

   interface
      subroutine regolith_push(toplayer,newlayer)
      use module_globals
      implicit none
      type(regolayertype),pointer :: toplayer
      type(regolayertype),intent(in) :: newlayer
      end subroutine regolith_push
   end interface

   interface
      subroutine regolith_pop(toplayer,oldlayer)
      use module_globals
      implicit none
      type(regolayertype),pointer :: toplayer
      type(regolayertype),intent(out) :: oldlayer
      end subroutine regolith_pop
   end interface

end module
