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
subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,totseb,turnover,dmix)
!subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,totseb,turnover)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_lineseg
   implicit none
   ! arguemnts
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(in) :: thetast,ri,rip1,zmin,zmax,erad,eradi,deltar
   type(regolayertype),intent(inout) :: newlayer
   real(DP),intent(inout) :: vmare,totseb
   logical,intent(inout) :: turnover
   real(DP),intent(inout) :: dmix 
   ! internal variables
   real(DP),parameter :: a = 0.936457 
   real(DP),parameter :: b = 1.12368
   type(regolayertype),pointer :: current
   real(DP) :: z,zstart,zend,rstart,rend,r
   real(DP) :: vsgly,x

   ! * Mixing
   real(DP) :: zmix 

   current => surfi%regolayer
   z = current%thickness
   zmix = z 
   zstart = 0.0_DP
   zend = z
 
   ! 1. Interacted 
   !if (zend >= zmax .and. zstart <= zmin) then
   !   vsgly = newlayer%thickness * user%pix**2
   !   vmare = vsgly * current%comp
   !   totseb = vsgly
      !write(*,*) 'l0',z,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,zmin,zmax,zstart,zend
   ! 2. Will be interacting 
   !else

     if (thetast>=0._DP) then
         rstart = ri
         rend = ri
     else 
         rstart = rip1
         rend = rip1
     end if

     do

       if (.not. associated(current%next)) exit !it should exit until it hit the very bottom.

       if (zend <= zmin) then
          if (zmax <= zend + current%next%thickness) then 
             vsgly = newlayer%thickness * user%pix**2
             vmare = vsgly * current%next%comp
             totseb = vsgly 
             !write(*,*) 'z<zmin',z,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,zmin,zmax,zstart,zend

             if (zend > zmix) then
                turnover = .true.
                dmix = dmix + vsgly / (user%pix * user%pix)
             end if

             exit
          else
             current => current%next
             z = z + current%thickness
             zstart = zend
             zend = z
          end if
       else if (zend > zmin .and. zend < zmax) then 
               rend = regolith_quadratic_func(zend,erad,ri,rip1,rstart)
               vsgly = 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b/eradi * rend) &
                       - tan(b/eradi * rstart)) - abs(b/eradi * rend - b/eradi * rstart))
               vmare = vmare + vsgly * current%comp
               totseb = totseb + vsgly

               if (zend > zmix) then
                  turnover = .true.
                  dmix = dmix + vsgly / (user%pix * user%pix)
               end if

               !write(*,*) 'lmid',z,current%comp,deltar,vsgly/user%pix**2,ri,rip1,rstart,rend,zmin,zmax,zstart,zend 
               current => current%next
               z = z + current%thickness
               r = rstart 
               rstart = rend 
               zstart = zend 
               zend = z
       else if (zend >= zmax .and. zstart <= zmin) then
               vsgly = newlayer%thickness * user%pix**2
               vmare = vsgly * current%comp
               totseb = vsgly
               !write(*,*) 'z>zmax',z,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,zmin,zmax,zstart,zend

               if (zstart > zmix) then
                  turnover = .true.
                  dmix = dmix + vsgly / (user%pix * user%pix)
               end if

               exit
       else if (zend >= zmax .and. zstart > zmin) then 
                ! last part of a stream tube
               vsgly = 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b/eradi * rip1) - tan(b/eradi * ri)) - & 
                       abs(b/eradi * rip1 - b/eradi * ri))
               vmare = vmare + (vsgly - totseb) * current%comp
               totseb = vsgly 

               if (zstart > zmix) then
                  turnover = .true.
                  dmix = dmix + vsgly / (user%pix * user%pix)
               end if

               !write(*,*) 'llast',z,current%comp,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,rstart,rend,zmin,zmax,zstart,zend
               exit
       !else 
       !        write(*,*) '??',ri,rip1,zmin,zmax,zstart,zend
       !        exit
       end if

      end do
      
   !end if

   !x = vmare/totseb
   !if (x /= x) write(*,*) 'linesegment summary',thetast,erad,ri,rip1,zstart,zend,zmin,zmax,vsgly,vmare,totseb
   !if (x > 1.1 .or. x < -1.d-15) then
   !   write(*,*) 'linseg, comp - x',x,erad,ri,rip1,zstart,zend,zmin,zmax,vsgly,vmare,totseb
   !   do i=1,50
   !      write(*,*) i,thick(i),r1(i),r2(i),z1(i),thick1(i)
   !   end do
   !end if

   return
end subroutine regolith_streamtube_lineseg
