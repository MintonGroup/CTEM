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
subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,vseg,newlayer,rm)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_traverse_streamtube
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   real(DP),intent(in)            :: deltar,ri,rip1,eradi,erado,vseg,rm
   type(regodatatype),intent(inout) :: newlayer

   ! Traversing a linked list 
   real(DP) :: zri,zrip1,cosi,coso,rzmax
   real(DP) :: erad,z,zmin,zmax,thetast

   erad = (eradi + erado)/2.0
   rzmax = erad * sqrt(3.0)/4.0
   cosi = regolith_quartic_func(ri,erad)
   zri  = erad * (1.0 - cosi) * cosi 
   coso = regolith_quartic_func(rip1,erad)
   zrip1 = erad * (1.0 - coso) * coso
   thetast = atan((zrip1-zri)/(rip1-ri))/PI*180.0

   zmin = min(zri,zrip1)
   zmax = max(zri,zrip1)

   !if (ri <= 0._DP .and. (rip1>rzmax .and. ri < rzmax)) then
   if (rip1>rzmax .and. ri < rzmax) then
      zmax = max(zmax,erad/4.0)
   end if

   z = surfi%regolayer%regodata%thickness

   if (z>=zmax) then 

      newlayer%thickness = vseg
      newlayer%comp      = vseg * surfi%regolayer%regodata%comp

   else 

     newlayer%thickness  = 0._DP
     newlayer%comp       = 0._DP
     call regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,vseg,newlayer,rm)

   end if

   return
end subroutine regolith_traverse_streamtube
