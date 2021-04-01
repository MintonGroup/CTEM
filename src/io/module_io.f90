!**********************************************************************************************************************************
!
!  Unit Name   : module_io
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for generic I/O 
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_io
use module_globals
implicit none
public 
save

   interface
      subroutine io_input(infile,user)
      use module_globals
      implicit none
      character(*), intent(in)  :: infile
      type(usertype),intent(inout) :: user
      end subroutine io_input
   end interface

   interface
      subroutine io_read_const(totalimpacts,ncount,curyear,restart,fracdone,masstot,seedarr)
      use module_globals
      implicit none
      integer(I8B),intent(out) :: totalimpacts
      integer(I4B),dimension(:),intent(out) :: seedarr
      integer(I4B),intent(out) :: ncount
      logical,intent(out) :: restart
      real(DP),intent(out) :: curyear,fracdone,masstot
      end subroutine io_read_const
   end interface

   interface
      subroutine io_write_const(totalimpacts,ncount,curyear,restart,fracdone,masstot,seedarr)
      use module_globals
      implicit none
      integer(I8B),intent(in) :: totalimpacts
      integer(I4B),intent(in) :: ncount
      logical,intent(in) :: restart
      real(DP),intent(in) :: curyear,fracdone,masstot
      integer(I4B),dimension(:),intent(in) :: seedarr
      end subroutine io_write_const
   end interface

   interface   
      subroutine io_read_surf(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(out) :: surf
      end subroutine io_read_surf
   end interface

  interface
      subroutine io_read_regotrack(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(out) :: surf
      end subroutine io_read_regotrack
  end interface  

  interface
      subroutine io_read_porotrack(user,surf)
		use module_globals
		implicit none
		type(usertype),intent(in) :: user
		type(surftype),dimension(:,:),intent(out) :: surf  
      end subroutine io_read_porotrack
  end interface  

  interface   
      subroutine io_write_surf(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf
      end subroutine io_write_surf
   end interface

  interface 
      subroutine io_write_regotrack(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf      
      end subroutine io_write_regotrack
  end interface 

  interface 
      subroutine io_write_porotrack(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf      
      end subroutine io_write_porotrack
  end interface 

  interface   
      subroutine io_crater_profile(user,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf
      end subroutine io_crater_profile
   end interface 

   interface
      subroutine io_read_prod(prod,user,domain)
      use module_globals
      implicit none
      real(DP),dimension(:,:),intent(out) :: prod
      type(usertype),intent(in) :: user
      type(domaintype),intent(in) :: domain
      end subroutine io_read_prod
   end interface

   interface
      subroutine io_read_craterlist(user,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(domaintype),intent(in) :: domain
      end subroutine io_read_craterlist
   end interface

   interface
      subroutine io_read_vdist(vdist,user,domain)
      use module_globals
      implicit none
      real(DP),dimension(:,:),intent(out) :: vdist
      type(usertype),intent(in) :: user
      type(domaintype),intent(in) :: domain
      end subroutine io_read_vdist
   end interface

   interface
      subroutine io_write_dist(pdist,crtscl,domain,mass)
      use module_globals
      implicit none
      real(DP),dimension(:,:),intent(in) :: pdist,crtscl
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: mass
      end subroutine io_write_dist
   end interface

   interface
      subroutine io_write_tally(tdist,tlist,odist,olist,oposlist,depthdiam)
      use module_globals
      implicit none
      real(DP),dimension(:,:),intent(in) :: tdist,tlist,odist
      real(DP),dimension(:),intent(in) :: olist
      real(SP),dimension(:,:),intent(in) :: oposlist
      real(SP),dimension(:),intent(in) :: depthdiam
      end subroutine io_write_tally
   end interface

   interface
      subroutine io_get_token(buffer, ilength, ifirst, ilast, ierr)
      use module_globals
      implicit none
      integer(I4B), intent(in)    :: ilength
      integer(I4B), intent(inout) :: ifirst
      integer(I4B), intent(out)   :: ilast, ierr
      character(*), intent(in)    :: buffer
      end subroutine io_get_token
   end interface

   interface
      subroutine io_ejecta_table(crater,domain,ejb,ejtble,filename)
      use module_globals
      implicit none
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      type(ejbtype),dimension(EJBTABSIZE),intent(in) :: ejb
      integer(I4B),intent(in) :: ejtble
      character(*),intent(in) :: filename
      end subroutine io_ejecta_table
   end interface

   interface
      subroutine io_updatePbar(message)
      use module_globals
      implicit none
      character(len=*),intent(in) :: message
      end subroutine io_updatePbar
   end interface

   interface
      subroutine io_resetPbar()
      use module_globals
      implicit none
      end subroutine io_resetPbar
   end interface

   interface
      subroutine io_splash()
      use module_globals
      implicit none
      end subroutine io_splash
   end interface

  !interface
  !    subroutine io_write_regodist(user,surf)
  !    use module_globals
  !    implicit none
  !    type(usertype),intent(in) :: user
  !    type(surftype),dimension(:,:),intent(in) :: surf
  !    end subroutine io_write_regodist
  !end interface

end module module_io
