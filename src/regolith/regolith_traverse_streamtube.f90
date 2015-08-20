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
subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,newlayer,vmare,totseb,turnover,dmix)
!subroutine regolith_traverse_streamtube(user,surfi,deltar,ri,rip1,eradi,erado,newlayer,vmare,totseb,turnover)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_traverse_streamtube
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   real(DP),intent(in)            :: deltar,ri,rip1,eradi,erado
   type(regolayertype),intent(inout) :: newlayer
   real(DP),intent(out)            :: vmare,totseb
   logical,intent(inout)           :: turnover
   real(DP),intent(inout)          :: dmix 
   !real(DP),dimension(200),intent(out) :: tots
   !real(DP),intent(out) :: thetast
   !integer(I4B),intent(out) :: cnt

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

   z = surfi%regolayer%thickness
   vmare = 0._DP
   totseb = 0._DP

   if (z>=zmax) then 
      vmare = newlayer%thickness * user%pix**2 * surfi%regolayer%comp
      totseb = newlayer%thickness * user%pix**2 
      !write(*,*) 'z>zmax',ri,rip1,erad/4.0,zmax,vmare/(user%pix**2),totseb/(user%pix**2)
   else 
     !if (ri == 0.0 .or. rip1>=eradi .or. abs(thetast)>10.0) then
     !call regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,totseb,&
     !     turnover)
     call regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,totseb,&
          turnover,dmix)
     !write(*,*) 'line',ri,rip1,zmax,vmare/(user%pix**2), totseb/(user%pix**2)
     !else
     !   call regolith_streamtube_cylinder(user,surfi,cosi,coso,ri,rip1,erad,eradi,deltar,thetast,vmare,totseb,&
     !        turnover,dmix)
        !write(*,*) 'cyli',ri,rip1,zmax,vmare/(user%pix**2), totseb/(user%pix**2)
     !end if
   end if

   return
end subroutine regolith_traverse_streamtube
