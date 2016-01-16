!**********************************************************************************************************************************
!
!  Unit Name   : module_seismic
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for seismic shaking 
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_seismic
use module_globals
implicit none
public 
save

   interface
      subroutine seismic_shake(user,surf,crater,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      end subroutine seismic_shake
   end interface

   interface
      subroutine seismic_distance(user,domain,crater,seisdis,maxhits,firstrun)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(domaintype),intent(in) :: domain
      type(cratertype),intent(in) :: crater
      real(DP),intent(out) :: seisdis
      integer(I4B),intent(out) :: maxhits
      logical,intent(inout) :: firstrun
      end subroutine seismic_distance
   end interface

   interface
      subroutine seismic_table_define(user,domain,crater,srim,totdiff,diffmax,srng,kdifflim,diffcnt,firstrun)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(domaintype),intent(in) :: domain
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: srim,totdiff,diffmax
      integer(I4B),intent(out) :: diffcnt
      integer(I4B),dimension(:),allocatable,intent(out) :: srng
      real(DP),dimension(:),allocatable,intent(out) :: kdifflim
      logical,intent(inout) :: firstrun
      end subroutine seismic_table_define
   end interface

   interface
      function seismic_kdiff_func(user,crater,lrad,gratio,firstrun,invflag) result(kdiff)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: lrad
      real(DP),intent(out) :: gratio
      logical,intent(inout) :: firstrun
      logical,intent(in) :: invflag
      real(DP) :: kdiff
      end function seismic_kdiff_func
   end interface

end module
