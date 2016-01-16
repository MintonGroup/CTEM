!**********************************************************************************************************************************
!
!  Unit Name   : module_init
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for initialization routines 
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_init
use module_globals
implicit none
public 
save

   interface
      subroutine init_surf(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(out) :: surf
      end subroutine init_surf
   end interface

   interface
      subroutine init_dist(user,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(domaintype),intent(inout) :: domain
      end subroutine init_dist
   end interface

   interface
      subroutine init_domain(user,crater,domain,prod,pdist,vdist,crtscl)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(inout) :: domain
      real(DP),dimension(:,:),intent(inout) :: prod,vdist
      real(DP),dimension(:,:),intent(out) :: pdist,crtscl
      end subroutine init_domain
   end interface 
  
   interface
      subroutine init_regolith_stack(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      end subroutine init_regolith_stack
   end interface

end module module_init
