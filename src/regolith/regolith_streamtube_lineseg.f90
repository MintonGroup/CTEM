!**********************************************************************************************************************************
!
!  Unit Name   : regolith_streamtube_lineseg
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : line segments approximation 
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
subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,totseb,&
           age_collector,xmints,xsfints,depthb,meltinejecta,totvol)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_lineseg
   implicit none
   ! arguemnts
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(in) :: thetast,ri,rip1,zmin,zmax,erad,eradi,deltar
   type(regodatatype),intent(inout) :: newlayer
   real(DP),intent(inout)           :: meltinejecta,totvol
   real(DP),intent(inout) :: vmare,totseb
   real(SP),dimension(:),intent(inout) :: age_collector
   real(DP),intent(in)             :: xmints
   real(DP),intent(in)             :: xsfints, depthb

   ! internal variables
   !type(regolisttype),pointer :: current
   type(regodatatype),dimension(:),allocatable :: current
   real(DP) :: z,zstart,zend,rstart,rend,r
   real(DP) :: vsgly,x,vseg, vsh
   integer(I4B) :: N

   ! Melt zone
   real(DP) :: recyratio, xsh, rst
   real(DP) :: theta1, theta2, r1, r2, vol
   real(DP) :: linmelt
   ! Shock damaged zone
   real(DP) :: ebh_recyl
 
   !current => surfi%regolayer
   allocate(current,source=surfi%regolayer)
   N = size(current)
   z = current(N)%thickness
   zstart = 0.0_DP
   zend = z
 
   if (thetast>=0._DP) then
      rstart = ri
      rend = ri
   else 
      rstart = rip1
      rend = rip1
   end if

   vol = 0.0_DP

   do

       !if (.not. associated(current%next)) exit !it should exit until it hit the very bottom.

       if (zend <= zmin) then

          if (zmax <= zend + current(N-1)%thickness) then 

             vsgly = newlayer%thickness * user%pix**2
             vmare = vsgly * current(N-1)%comp
             totseb = vsgly 

             if (ri > xmints) then
                vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,ri,rip1)
                ebh_recyl        = (1.0 - newlayer%meltfrac) * newlayer%thickness - vsh / user%pix**2
                recyratio        = max(ebh_recyl,0.0_DP) / current(N-1)%thickness
                age_collector(:) = age_collector(:) + current(N-1)%age(:) * recyratio
                vol              = vol + sum(current(N-1)%age(:)) * recyratio
                linmelt = current(N-1)%meltfrac * vsgly * recyratio
                meltinejecta = meltinejecta + linmelt
                totvol = totvol + vsgly
             else if (ri <= xmints .and. rip1 > xmints) then 
                     vseg             = regolith_streamtube_volume_func(eradi,xmints,rip1,deltar)
                     vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,ri,rip1) 
                     recyratio        = max((vseg-vsh),0.0_DP) / (user%pix**2) / current(N-1)%thickness 
                     age_collector(:) = age_collector(:) + current(N-1)%age(:) * recyratio
                     vol              = vol + sum(current(N-1)%age(:)) * recyratio
                     linmelt = current(N-1)%meltfrac * vseg * recyratio
                     meltinejecta = meltinejecta + linmelt
                     totvol = totvol + vseg
             end if
             exit
          else

             !current => current%next
            N = N - 1
            z = z + current(N)%thickness
            zstart = zend
            zend = z

          end if

       else if (zend > zmin .and. zend < zmax) then 

               rend   = regolith_quadratic_func(zend,erad,ri,rip1,rstart)
               vsgly  = regolith_streamtube_volume_func(eradi,rstart,rend,deltar) 
               vmare  = vmare + vsgly * current(N)%comp
               totseb = totseb + vsgly

               ! A segment coming from the side of impact site, ri
               if (thetast>=0_DP .and. rend > xmints) then 
                  vseg             = regolith_streamtube_volume_func(eradi,max(xmints,rstart),rend,deltar)
                  vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,rstart,rend)
                  recyratio        = max((vseg-vsh),0.0_DP) / (user%pix**2) / current(N)%thickness
                  age_collector(:) = age_collector(:) + current(N)%age(:) * recyratio
                  vol              = vol + sum(current(N)%age(:)) * recyratio
                  linmelt = current(N)%meltfrac * vseg * recyratio
                  meltinejecta = meltinejecta + linmelt
                  totvol = totvol + vseg
               end if

               ! A segment coming from the side of emerging location of a streamtube rip1 
               if (thetast<0._DP .and. rstart > xmints) then
                  vseg             = regolith_streamtube_volume_func(eradi,max(xmints,rend),rstart,deltar)
                  vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,rend,rstart)
                  recyratio        = max((vseg-vsh),0.0_DP) / (user%pix**2) / current(N)%thickness
                  age_collector(:) = age_collector(:) + current(N)%age(:) * recyratio
                  vol              = vol + sum(current(N)%age(:)) * recyratio
                  linmelt = current(N)%meltfrac * vseg * recyratio
                  meltinejecta = meltinejecta + linmelt
                  totvol = totvol + vseg
               end if
               !current => current%next
               N = N - 1
               z = z + current(N)%thickness
               r = rstart 
               rstart = rend 
               zstart = zend 
               zend = z

       else if (zend >= zmax .and. zstart <= zmin) then

               vsgly = newlayer%thickness * user%pix**2
               vmare = vsgly * current(N)%comp
               totseb = vsgly

               if (rip1 > xmints) then
                  vseg             = regolith_streamtube_volume_func(eradi,max(xmints,ri),rip1,deltar)
                  vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,ri,rip1)
                  recyratio        = max((vseg-vsh),0.0_DP) / (user%pix**2) / current(N)%thickness
                  age_collector(:) = age_collector(:) + current(N)%age(:) * recyratio
                  vol              = vol + sum(current(N)%age(:)) * recyratio
                  linmelt = current(N)%meltfrac * vseg * recyratio
                  meltinejecta = meltinejecta + linmelt
                  totvol = totvol + vseg
               end if

               exit

       else if (zend >= zmax .and. zstart > zmin) then 
               ! last part of a stream tube
               vsgly = regolith_streamtube_volume_func(eradi,ri,rip1,deltar)
               vmare = vmare + (vsgly - totseb) * current(N)%comp
               totseb = vsgly
               linmelt = current(N)%meltfrac * vsgly
               meltinejecta = meltinejecta + linmelt
               totvol = totvol + vsgly
               exit
       end if

   end do
      
   return
end subroutine regolith_streamtube_lineseg
