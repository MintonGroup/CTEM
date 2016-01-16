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
      subroutine regolith_push(surf,newlayer)
      use module_globals
      implicit none
      type(surftype),intent(inout) :: surf
      type(regolayertype),intent(in) :: newlayer
      end subroutine regolith_push
   end interface

   interface
      subroutine regolith_pop(surfi)
      use module_globals
      implicit none
      type(surftype),intent(inout):: surfi
      end subroutine regolith_pop
   end interface

   interface
      subroutine regolith_traverse_pop(elchange,surfi)
      use module_globals
      implicit none
      real(DP),intent(in)         :: elchange
      type(surftype),intent(inout):: surfi      
      end subroutine 
   end interface

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
      subroutine regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,comp)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(in) :: ejtble
      type(ejbtype),dimension(ejtble),intent(in)   :: ejb
      real(DP),intent(in)          :: xp,yp,lrad,ebh
      real(DP),intent(out)         :: comp 
      integer(I4B),intent(in)      :: xpi,ypi
      end subroutine regolith_streamtube
   end interface

   interface 
      subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,newlayer,vmare,totseb,turnover,dmix)
!      subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,newlayer,vmare,totseb,turnover)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      real(DP),intent(in)            :: deltar,ri,rip1,eradi,erado
      type(regolayertype),intent(inout) :: newlayer
      real(DP),intent(out) :: vmare,totseb
      logical,intent(inout) :: turnover
      real(DP),intent(inout) :: dmix 
      end subroutine regolith_traverse_streamtube
   end interface

   interface 
      subroutine regolith_subpixel_streamtube(user,surfi,deltar,ri,rip1,eradi,newlayer,vmare,totseb,turnover,dmix)
!      subroutine regolith_subpixel_streamtube(user,surfi,deltar,ri,rip1,eradi,newlayer,vmare,totseb,turnover)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(inout) :: surfi
      real(DP),intent(in)            :: deltar,ri,rip1,eradi
      type(regolayertype),intent(inout) :: newlayer
      real(DP),intent(out) :: vmare,totseb
      logical, intent(inout) :: turnover
      real(DP),intent(inout) :: dmix
      end subroutine regolith_subpixel_streamtube
   end interface

   interface 
      subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,&
      totseb,turnover,dmix)
!      subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,&
!      totseb,turnover)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(in) :: surfi
      real(DP),intent(in) :: thetast,ri,rip1,zmin,zmax,erad,eradi,deltar
      type(regolayertype),intent(inout) :: newlayer
      real(DP),intent(inout) :: vmare,totseb
      logical,intent(inout) :: turnover
      real(DP),intent(inout) :: dmix 
      end subroutine regolith_streamtube_lineseg
   end interface 

!   interface 
!      subroutine regolith_streamtube_cylinder(user,surfi,cosi,coso,ri,rip1,erad,eradi,deltar,thetast,vmare,totseb,turnover,dmix)
!      subroutine regolith_streamtube_cylinder(user,surfi,cosi,coso,ri,rip1,erad,eradi,deltar,thetast,vmare,totseb,turnover)
!      use module_globals 
!      implicit none
!      type(usertype),intent(in) :: user
!      type(surftype),intent(inout) :: surfi
!      real(DP),intent(in) :: cosi,coso,ri,rip1,erad,eradi,deltar,thetast
!      real(DP),intent(inout) :: vmare,totseb
!      logical,intent(inout) :: turnover
!      real(DP),intent(inout) :: dmix 
!      end subroutine regolith_streamtube_cylinder
!   end interface

   interface 
      subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots,turnover,dmix)
!      subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots,turnover)
      use module_globals 
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),intent(in) :: surfi
      real(DP),intent(in) :: deltar
      real(DP),intent(inout) :: totmare,tots
      logical,intent(inout)  :: turnover
      real(DP),intent(inout) :: dmix 
      end subroutine regolith_streamtube_head
   end interface

   interface
      subroutine regolith_rays(user,crater,domain,ejtble,ejb)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(cratertype),intent(in) :: crater
      type(domaintype),intent(in) :: domain
      integer(I4B),intent(in) :: ejtble
      type(ejbtype),dimension(ejtble),intent(in)   :: ejb
      end subroutine regolith_rays
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
      subroutine regolith_monte_carlo_layer(surfi,erad,deltar,mceb,compmc)
      use module_globals 
      implicit none
      type(surftype),intent(in)    :: surfi
      real(DP),intent(in)          :: erad,deltar
      real(DP),intent(out)         :: mceb,compmc
      end subroutine regolith_monte_carlo_layer
   end interface

   interface 
      subroutine regolith_streamtube_mc(surfi,ri,rip1,erad,deltar,totmc,cnt)
      use module_globals 
      implicit none
      type(surftype),intent(in)      :: surfi
      real(DP),intent(in)            :: ri,rip1,erad,deltar
!      real(DP),intent(inout) :: totmc
      real(DP),dimension(100),intent(out) :: totmc
      integer,intent(out)            :: cnt
      end subroutine regolith_streamtube_mc
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
      subroutine regolith_reworking_zone(user,surf,finterval)
      use module_globals
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      real(DP),intent(in) :: finterval
      end subroutine regolith_reworking_zone
   end interface

end module
