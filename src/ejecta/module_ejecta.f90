!****h* ejecta/module_ejecta
! Name
!   module_ejecta -- Module for ejecta
! Notes
!   It includes crater ray model based on Superformula. 
!***

!**********************************************************************************************************************************
!
!  Unit Name   : module_crater
!  Unit Type   : module
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for crater creation, including scaling laws and Monte Carlo
!                routines
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_ejecta
use module_globals
implicit none
public 
save

   interface
      subroutine ejecta_emplace(user,surf,crater,domain,ejb,ejtble)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(in) :: ejtble
      type(ejbtype),dimension(ejtble),intent(in)   :: ejb
      end subroutine ejecta_emplace
   end interface

   interface
      subroutine ejecta_rootfind(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(inout) :: erad
      real(DP),intent(in) :: lrad
      real(DP),intent(out) :: vejsq,ejang
      logical,intent(inout) :: firstrun
      end subroutine ejecta_rootfind
   end interface

   interface
      subroutine ejecta_blanket(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: erad
      real(DP),intent(out) :: lrad,vejsq,ejang
      logical,intent(inout) :: firstrun
      end subroutine ejecta_blanket
   end interface

   interface
      subroutine ejecta_thickness(user,crater,erad1,erad2,lrad1,lrad2,thick)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: erad1,erad2,lrad1,lrad2
      real(DP),intent(out) :: thick
      end subroutine ejecta_thickness
   end interface

   interface
      function ejecta_blanket_func(user,crater,domain,erad,lrad,vejsq,ejang,firstrun) result(ans)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: erad,lrad
      logical,intent(inout) :: firstrun
      real(DP),intent(out) :: vejsq,ejang
      real(DP) :: ans
      end function ejecta_blanket_func
   end interface

   interface
      subroutine ejecta_table_define(user,crater,domain,ejb,ejtble,melt)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(inout) :: domain
      type(ejbtype),dimension(EJBTABSIZE),intent(out) :: ejb
      integer(I4B),intent(out) :: ejtble
      real(DP),intent(out),optional :: melt
      end subroutine ejecta_table_define
   end interface

   interface
      subroutine ejecta_interpolate(crater,domain,lrad,ejb,ejtble,ebh,vsq,theta,melt)
      use module_globals
      implicit none
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in)  :: lrad
      integer(I4B),intent(in) :: ejtble
      type(ejbtype),dimension(ejtble),intent(in) :: ejb
      real(DP),intent(out) :: ebh
      real(DP),intent(out),optional :: vsq,theta
      real(DP),intent(out),optional :: melt
      end subroutine ejecta_interpolate
   end interface

   interface
      subroutine ejecta_soften(user,surf,N,indarray,cumulative_elchange)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      integer(I4B),intent(in) :: N
      integer(I4B),dimension(2,N,N),intent(in) :: indarray
      real(DP),dimension(N,N),intent(inout) :: cumulative_elchange 
      end subroutine ejecta_soften
   end interface

   interface
      subroutine ejecta_distance_estimate(user,crater,domain,ejdis_estimate)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in)    :: domain
      real(DP),intent(out) :: ejdis_estimate
      end subroutine ejecta_distance_estimate
   end interface
end module
