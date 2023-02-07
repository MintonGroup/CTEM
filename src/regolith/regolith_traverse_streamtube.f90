!**********************************************************************************************************************************
!
!  Unit Name   : regolith_streamtube_volume
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Summing up a history of a stream tube 
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,newlayer,vmare,totseb,&
           age_collector,xmints,xsfints,depthb,meltinejecta,totvol)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_traverse_streamtube
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   real(DP),intent(in)            :: deltar,ri,rip1,eradi,erado
   type(regodatatype),intent(inout) :: newlayer
   real(DP),intent(inout)            :: meltinejecta,totvol
   real(DP),intent(out)            :: vmare,totseb
   real(SP),dimension(:),intent(inout) :: age_collector
   real(DP),intent(in)             :: xmints
   real(DP),intent(in)             :: xsfints, depthb

   ! Internal variables

   ! Traversing a linked list 
   real(DP) :: zri,zrip1,cosi,coso,rzmax
   real(DP) :: erad,z,zmin,zmax,thetast,vseg
   integer(I4B) :: N

   real(DP),parameter :: a = 0.936457
   real(DP),parameter :: b = 1.12368

   ! Melt zone 
   real(DP) :: recyratio

   ! Shock zone
   real(DP) :: vsh

   !executable code

   erad = (eradi + erado)/2.0
   rzmax = erad * sqrt(3.0)/4.0
   if (ri .eq. 0) then
      cosi = 0.0_DP
   else
      cosi = regolith_quartic_func(ri,erad)
   end if
   zri  = erad * (1.0 - cosi) * cosi 
   coso = regolith_quartic_func(rip1,erad)
   zrip1 = erad * (1.0 - coso) * coso
   thetast = atan2((zrip1-zri),(rip1-ri))/PI*180.0

   zmin = min(zri,zrip1)
   zmax = max(zri,zrip1)

   !if (ri <= 0._DP .and. (rip1>rzmax .and. ri < rzmax)) then
   if (rip1>rzmax .and. ri < rzmax) then
      zmax = max(zmax,erad/4.0)
   end if

   N = size(surfi%regolayer)

   z = surfi%regolayer(N)%thickness
   vmare  = 0._DP
   totseb = 0._DP

   if (z>=zmax) then 

      vmare = newlayer%thickness * user%pix**2 * surfi%regolayer(N)%comp
      totseb = newlayer%thickness * user%pix**2 
      if (rip1 > xmints .and. ri < xmints) then
         vseg             = regolith_streamtube_volume_func(eradi,max(xmints,ri),rip1,deltar)
         vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,ri,rip1)
         recyratio        = max(vseg-vsh,0.0_DP) / (user%pix**2) / surfi%regolayer(N)%thickness
         age_collector(:) = age_collector(:) + surfi%regolayer(N)%age(:) * recyratio
         meltinejecta     = meltinejecta + surfi%regolayer(N)%meltfrac * vseg * recyratio
         totvol = totvol + vseg
      end if

   else 

     call regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,&
          newlayer,vmare,totseb,age_collector,xmints,xsfints,depthb,meltinejecta,totvol)

   end if

   return
end subroutine regolith_traverse_streamtube
