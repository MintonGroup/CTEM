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
      subroutine regolith_melt_zone(user,crater,dimp,vimp,rmelt,depthb)
      use module_globals
      type(usertype),intent(in) :: user
      type(cratertype),intent(inout) :: crater
      real(DP),intent(in)  :: dimp,vimp
      real(DP),intent(out) :: rmelt,depthb
      end subroutine regolith_melt_zone
   end interface

   interface 
      subroutine regolith_melt_fraction(dimp,depthb,erad1,erad2,rmelt,meltfrac)
      use module_globals
      real(DP),intent(in)       :: dimp,depthb,erad1,erad2,rmelt
      real(DP),intent(inout)    :: meltfrac
      end subroutine regolith_melt_fraction
   end interface

   interface 
      function regolith_melt_func(r,depthb,erad) result(vol)
      use module_globals
      real(DP),intent(in) :: r,depthb,erad 
      real(DP) :: vol
      end function regolith_melt_func
   end interface

   interface 
      subroutine regolith_transport(user,surfi,crater,domain,ejb,ejtble,lrad,ebh,comp)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(in) :: ejtble
      type(ejbtype),dimension(ejtble),intent(in)   :: ejb
      real(DP),intent(in)          :: lrad,ebh,comp
      end subroutine regolith_transport
   end interface

   interface 
      subroutine regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,comp,eradc)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(in) :: ejtble
      type(ejbtype),dimension(ejtble),intent(in)   :: ejb
      real(DP),intent(in)          :: xp,yp,lrad,ebh
      real(DP),intent(out)         :: comp, eradc
      integer(I4B),intent(in)      :: xpi,ypi
      end subroutine regolith_streamtube
   end interface

   interface 
      subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,newlayer,vmare,totseb)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      real(DP),intent(in)            :: deltar,ri,rip1,eradi,erado
      type(regodatatype),intent(inout) :: newlayer
      real(DP),intent(out) :: vmare,totseb
      end subroutine regolith_traverse_streamtube
   end interface

   interface 
      subroutine regolith_subpixel_streamtube(user,surfi,deltar,ri,rip1,eradi,newlayer,vmare,totseb)      
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      real(DP),intent(in)            :: deltar,ri,rip1,eradi
      type(regodatatype),intent(inout) :: newlayer
      real(DP),intent(out) :: vmare,totseb
      end subroutine regolith_subpixel_streamtube
   end interface

   interface 
      subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,&
      totseb)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(in) :: surfi
      real(DP),intent(in) :: thetast,ri,rip1,zmin,zmax,erad,eradi,deltar
      type(regodatatype),intent(inout) :: newlayer
      real(DP),intent(inout) :: vmare,totseb
      end subroutine regolith_streamtube_lineseg
   end interface 

   interface 
      subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(in) :: surfi
      real(DP),intent(in) :: deltar
      real(DP),intent(inout) :: totmare,tots
      end subroutine regolith_streamtube_head
   end interface

   interface
      function regolith_circle_sector_func(deltar,zstart,zend) result(vhead)
      use module_globals 
      implicit none
      real(DP),intent(in) :: deltar,zstart,zend
      real(DP)            :: vhead
      end function regolith_circle_sector_func
   end interface

   interface 
      function regolith_quartic_func(rpj,r) result(z)
      use module_globals
      implicit none
      real(DP),intent(in)              :: rpj,r
      real(DP)                         :: z
      end function regolith_quartic_func
   end interface

   interface 
      function regolith_quadratic_func(z,erad,ri,rip1,rstart) result(r)
      use module_globals
      implicit none
      real(DP),intent(in)              :: z,erad,ri,rip1,rstart
      real(DP)          :: r
      end function regolith_quadratic_func
   end interface

   interface 
      function regolith_cylinder_func(z,zc,rcyl) result(aseg)
      use module_globals
      implicit none
      real(DP),intent(in)              :: z,zc,rcyl
      real(DP)                         :: aseg
      end function regolith_cylinder_func
   end interface

   interface 
      function regolith_cubic_func(c1,c2) result(deltar)
      use module_globals
      implicit none
      real(DP),intent(in)   :: c1,c2
      real(DP) :: deltar
      end function regolith_cubic_func
   end interface

   interface 
      subroutine regolith_subcrater_diffusion(user,surf,grad,mi,mj,mi2,mj2,regopop,regopush,comp,&
      n,roundoff)
      use module_globals 
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(in) :: surf
      real(DP),dimension(4),intent(in) :: grad
      integer(I4B),intent(in) :: mi,mj
      integer(I4B),dimension(4),intent(in) :: mi2,mj2
      real(DP),intent(inout) :: regopop,regopush,comp
      integer(I4B),intent(inout) :: n
      real(DP),intent(inout) :: roundoff
      end subroutine regolith_subcrater_diffusion
   end interface

   interface
      subroutine regolith_subcrater_mix(user,surf,domain,nflux,finterval,p)
      use module_globals
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(domaintype),intent(in) :: domain
      real(DP),dimension(:,:),intent(in) :: nflux ! impact rate (number of craters per m^2 per year)
      real(DP),intent(in) :: finterval  ! time elapsed ratio to the total time 
      real(DP),dimension(:,:),intent(in) :: p
      end subroutine regolith_subcrater_mix
   end interface

   interface
      subroutine regolith_mix(surf,d)
      use module_globals
      type(surftype),intent(inout) :: surf
      real(DP),intent(in) :: d
      end subroutine regolith_mix
   end interface

   interface
      subroutine regolith_depth_model(user,domain,finterval,nflux,p)
      use module_globals
      type(usertype),intent(in) :: user
      type(domaintype),intent(in) :: domain
      real(DP),intent(in) :: finterval
      real(DP),dimension(:,:),intent(in) :: nflux
      real(DP),dimension(:,:),intent(out) :: p
      end subroutine regolith_depth_model
   end interface

end module
