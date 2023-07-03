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
                                 mass,fracdone,nflux,ntotcrat,curyear,rclist)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout)  :: surf
      type(cratertype),intent(inout)               :: crater
      type(domaintype),intent(inout)               :: domain
      real(DP),dimension(:,:),intent(in)           :: prod,vdist
      integer(I8B),dimension(:),intent(inout)      :: production_list            
      integer(I4B),intent(out)                     :: ntrue
      integer(I4B),intent(out)                     :: vistrue
      integer(I4B),intent(out)                     :: ntotkilled
      real(DP),dimension(:,:),allocatable,intent(inout) :: truelist
      real(DP),intent(out)                         :: mass
      real(DP),intent(out)                         :: fracdone
      real(DP),dimension(:,:),intent(in)           :: nflux 
      integer(I8B),intent(in)                      :: ntotcrat
      real(DP),intent(in)                          :: curyear
      real(DP),dimension(:,:), intent(in)              :: rclist !array of 'real' craters for quasiMC
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
      subroutine crater_dimensions(user,crater,domain)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in)    :: domain
      end subroutine crater_dimensions
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
      subroutine crater_emplace(user,surf,crater,domain,deltaMtot,incval,nmeltsheet)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(inout) :: domain
      real(DP),intent(out) :: deltaMtot
      integer(I4B),intent(out) :: incval,nmeltsheet
      end subroutine crater_emplace
   end interface

   interface
      subroutine crater_realistic_topography(user,surf,crater,domain,ejecta_dem)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),dimension(:,:),intent(inout) :: ejecta_dem
      end subroutine crater_realistic_topography
   end interface

   interface
      subroutine crater_form_interior(user,surfi,crater,x_relative, y_relative, newelev,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(inout) :: crater
      real(DP),intent(in) :: x_relative, y_relative 
      real(DP),intent(in) :: newelev
      real(DP),intent(out) :: deltaMi
      end subroutine crater_form_interior
   end interface

   interface
      subroutine crater_form_exterior(user,surfi,crater,domain,lradsq,newelev,rimheight,deltaMi)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: newelev,lradsq,rimheight
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
      real(DP),dimension(:,:),intent(inout)  :: truelist
      real(DP),dimension(:,:),intent(out) :: truedist
      end subroutine crater_tally_true
   end interface

   interface
      subroutine crater_tally_observed(user,surf,domain,nkilled,onum,obsdist,obslist,oposlist,depthdiam,degradation_state)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(out) :: nkilled,onum
      real(DP),dimension(:,:),intent(out),optional  :: obsdist
      real(DP),dimension(:),intent(out),allocatable,optional :: obslist
      real(SP),dimension(:,:),intent(out),allocatable,optional :: oposlist
      real(SP),dimension(:),intent(out),allocatable,optional :: depthdiam
      real(DP),dimension(:),intent(out),allocatable,optional :: degradation_state
      end subroutine crater_tally_observed
   end interface

   interface
      subroutine crater_slope_collapse(user,surf,crater,domain,critical_value,deltaMtot)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: critical_value
      real(DP),intent(inout) :: deltaMtot
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
      subroutine crater_soften_accumulate(user,surf,crater,domain,kdiff)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      real(DP),dimension(:,:),intent(inout) :: kdiff
      end subroutine crater_soften_accumulate
   end interface
   interface
      subroutine crater_subpixel_diffusion(user,surf,nflux,domain,finterval,kdiffin)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      real(DP),dimension(:,:),intent(in) :: nflux
      type(domaintype),intent(in)    :: domain
      real(DP),intent(in) :: finterval
      real(DP),dimension(:,:),intent(inout) :: kdiffin
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
      function crater_degradation_function(user,r) result(Kd)
      use module_globals
      type(usertype),intent(in) :: user
      real(DP),intent(in) :: r
      real(DP) :: Kd
      end function crater_degradation_function
   end interface

   interface
      function crater_get_degradation_state(user,surf,crater,dd) result(Kval)
      use module_globals
      implicit none
      ! Arguments
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf
      type(cratertype),intent(in) :: crater
      real(SP),intent(out) :: dd
      ! Return value
      real(DP) :: Kval
      end function crater_get_degradation_state
   end interface

   interface
      function crater_visibility(user,crater,Kval) result(iscountable)
      use module_globals
      implicit none
      ! Arguments
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: Kval
      ! Result variable
      logical :: iscountable       
      end function crater_visibility
   end interface


   interface
      function crater_profile(user,crater,r) result(h)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      real(DP),intent(in) :: r
      real(DP) :: h
      end function crater_profile
   end interface


   interface
      function crater_profile_find_r_inner_wall(user,crater) result(r_inner_wall)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      real(DP) :: r_inner_wall
      end function crater_profile_find_r_inner_wall
   end interface

   interface
      subroutine crater_superdomain(user,surf,age,age_resolution,prod,nflux,domain,finterval)
      use module_globals
      type(usertype),intent(in)                           :: user
      type(surftype),dimension(:,:),intent(inout)         :: surf
      real(DP),intent(in)                                 :: age
      real(DP),intent(in)                                 :: age_resolution
      real(DP),dimension(:,:),intent(in)                  :: prod,nflux
      type(domaintype),intent(in)                         :: domain
      real(DP),intent(in)                                 :: finterval
      end subroutine crater_superdomain
   end interface

end module
