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
module module_crater
use module_globals
implicit none
public 
save

interface
   subroutine crater_mass_conservation(user,surf,crater)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(in)  :: crater
   end subroutine crater_mass_conservation
end interface


   interface
      subroutine crater_populate(user,surf,crater,domain,prod,production_list,vdist,ntrue,vistrue,ntotkilled,truelist,&
                                 mass,fracdone,nflux,ntotcrat)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout)  :: surf
      type(cratertype),intent(inout)               :: crater
      type(domaintype),intent(inout)               :: domain
      real(DP),dimension(:,:),intent(in)           :: prod,vdist
      integer(I8B),dimension(:),intent(inout)         :: production_list            
      integer(I4B),intent(out)                     :: ntrue
      integer(I4B),intent(out)                     :: vistrue
      integer(I4B),intent(out)                     :: ntotkilled
      real(DP),dimension(:,:),intent(out)          :: truelist
      real(DP),intent(out)                         :: mass
      real(DP),intent(out)                         :: fracdone
      real(DP),dimension(:,:),intent(in)           :: nflux 
      integer(I8B),intent(in)                      :: ntotcrat
      end subroutine crater_populate
   end interface

   interface
      subroutine crater_scale(user,psize,crad,grad,strflag,sinimpang,impvel)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      real(DP),intent(in) :: psize,sinimpang,impvel
      real(DP),intent(out) :: crad,grad
      integer(I4B),intent(in) :: strflag
      end subroutine crater_scale
   end interface

   interface
      subroutine crater_generate(user,crater,domain,prod,production_list,vdist,surf)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in)    :: domain
      real(DP),dimension(:,:),intent(in),optional :: prod,vdist
      integer(I8B),dimension(:),intent(inout),optional  :: production_list            
      type(surftype),dimension(:,:),intent(in),optional :: surf
      end subroutine crater_generate
   end interface

   interface
      subroutine crater_find_visible(user,crater,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in)    :: domain
      end subroutine crater_find_visible
   end interface

   interface
      subroutine crater_averages(user,surf,crater)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf
      type(cratertype),intent(inout) :: crater
      end subroutine crater_averages
   end interface

   interface
      subroutine crater_emplace(user,surf,crater,domain,ejbmass)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(inout) :: domain
      real(DP),intent(in) :: ejbmass
      end subroutine crater_emplace
   end interface

   interface
      subroutine crater_form_interior(user,surfi,crater,lradsq,newelev,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: lradsq
      real(DP),intent(in) :: newelev
      real(DP),intent(out) :: deltaMi
      end subroutine crater_form_interior
   end interface

   interface
      subroutine crater_form_exterior(user,surfi,crater,domain,lradsq,newelev,rd,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: newelev,lradsq,rd
      real(DP),intent(out) :: deltaMi
      end subroutine crater_form_exterior
   end interface

   interface
      function crater_form_exterior_func(user,surf,crater,domain,rd,deltaMtot,lastloop) result(ans)
      use module_globals
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: rd,deltaMtot
      logical,intent(in) :: lastloop
      real(DP) :: ans
      end function crater_form_exterior_func
   end interface

   interface
      subroutine crater_form_exterior_rootfind(user,surf,crater,domain,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) ::   deltaMtot
      end subroutine crater_form_exterior_rootfind
   end interface

   interface
      subroutine crater_record(user,surf,crater)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      end subroutine crater_record
   end interface

   interface
      subroutine crater_tally_true(domain,truelist,ntrue,truedist)
      use module_globals
      implicit none
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(in)   :: ntrue
      real(DP),dimension(TRUECOLS,ntrue),intent(inout)  :: truelist
      real(DP),dimension(:,:),intent(out) :: truedist
      end subroutine crater_tally_true
   end interface

   interface
      subroutine crater_tally_observed(user,surf,domain,nkilled,onum,obsdist,obslist,oposlist,&
                                       original_depth,current_depth,deviation_sigma,p)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(out) :: nkilled,onum
      real(DP),dimension(:,:),intent(out),optional  :: obsdist
      real(DP),dimension(:),intent(out),allocatable,optional :: obslist
      real(SP),dimension(:,:),intent(out),allocatable,optional :: oposlist
      real(SP),dimension(:),intent(out),allocatable,optional :: original_depth,current_depth,deviation_sigma,p
      end subroutine crater_tally_observed
   end interface

interface
   subroutine crater_tally_calibrated_count(user,diameter,current_depth,original_depth,&
                                            deviation_sigma,countable,killable,p)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   real(DP),intent(in) :: diameter
   real(SP),intent(in) :: current_depth,original_depth,deviation_sigma
   logical,intent(out) :: countable,killable
   real(SP),intent(out) :: p
   end subroutine crater_tally_calibrated_count
end interface

   interface
      subroutine crater_slope_collapse(user,surf,crater,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      end subroutine crater_slope_collapse
   end interface

   interface
      subroutine crater_soften(user,surf,crater,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      end subroutine crater_soften
   end interface


   interface
      subroutine crater_subpixel_diffusion(user,surf,prod,crtscl,domain,finterval)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      real(DP),dimension(:,:),intent(in) :: prod,crtscl
      type(domaintype),intent(in)    :: domain
      real(DP),intent(in) :: finterval
      end subroutine crater_subpixel_diffusion
   end interface

   interface
      subroutine crater_make_list(domain,prod,ntotcrat,production_list)
      use module_globals
      implicit none
      type(domaintype),intent(in) :: domain
      real(DP),dimension(:,:),intent(in) :: prod
      integer(I8B),intent(out) :: ntotcrat
      integer(I8B),dimension(:),intent(out) :: production_list 
      end subroutine crater_make_list
   end interface

   interface
      function crater_critical_slope(user,crater,iradsq) result(critical)
      use module_globals
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      integer(I4B),intent(in) :: iradsq
      real(DP) :: critical
      end function crater_critical_slope
   end interface

end module
