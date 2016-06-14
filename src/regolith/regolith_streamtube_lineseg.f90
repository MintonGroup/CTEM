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
subroutine regolith_streamtube_lineseg(user,surfi,thetast,ri,rip1,zmin,zmax,erad,eradi,deltar,newlayer,vmare,totseb)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_lineseg
   implicit none
   ! arguemnts
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(in) :: thetast,ri,rip1,zmin,zmax,erad,eradi,deltar
   type(regodatatype),intent(inout) :: newlayer
   real(DP),intent(inout) :: vmare,totseb
   ! internal variables
   real(DP),parameter :: a = 0.936457 
   real(DP),parameter :: b = 1.12368
   type(regolisttype),pointer :: current
   real(DP) :: z,zstart,zend,rstart,rend,r
   real(DP) :: vsgly,x

   ! * Mixing
   real(DP) :: zmix 

   current => surfi%regolayer
   z = current%regodata%thickness
   zmix = z 
   zstart = 0.0_DP
   zend = z
 
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
          if (zmax <= zend + current%next%regodata%thickness) then 
             vsgly = newlayer%thickness * user%pix**2
             vmare = vsgly * current%next%regodata%comp
             totseb = vsgly 
             !write(*,*) 'z<zmin',z,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,zmin,zmax,zstart,zend
             exit
          else
             current => current%next
             z = z + current%regodata%thickness
             zstart = zend
             zend = z
          end if
       else if (zend > zmin .and. zend < zmax) then 
               rend = regolith_quadratic_func(zend,erad,ri,rip1,rstart)
               vsgly = 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b/eradi * rend) &
                       - tan(b/eradi * rstart)) - abs(b/eradi * rend - b/eradi * rstart))
               vmare = vmare + vsgly * current%regodata%comp
               totseb = totseb + vsgly
               !write(*,*) 'lmid',z,current%comp,deltar,vsgly/user%pix**2,ri,rip1,rstart,rend,zmin,zmax,zstart,zend 
               current => current%next
               z = z + current%regodata%thickness
               r = rstart 
               rstart = rend 
               zstart = zend 
               zend = z
       else if (zend >= zmax .and. zstart <= zmin) then
               vsgly = newlayer%thickness * user%pix**2
               vmare = vsgly * current%regodata%comp
               totseb = vsgly
               !write(*,*) 'z>zmax',z,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,zmin,zmax,zstart,zend
               exit
       else if (zend >= zmax .and. zstart > zmin) then 
                ! last part of a stream tube
               vsgly = 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b/eradi * rip1) - tan(b/eradi * ri)) - & 
                       abs(b/eradi * rip1 - b/eradi * ri))
               vmare = vmare + (vsgly - totseb) * current%regodata%comp
               totseb = vsgly 
               !write(*,*) 'llast',z,current%comp,vmare/user%pix**2,vsgly/user%pix**2,ri,rip1,rstart,rend,zmin,zmax,zstart,zend
               exit
       !else 
       !        write(*,*) '??',ri,rip1,zmin,zmax,zstart,zend
       !        exit
       end if

   end do
      
   return
end subroutine regolith_streamtube_lineseg
