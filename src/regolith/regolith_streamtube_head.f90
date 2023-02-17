!**********************************************************************************************************************************
!
!  Unit Name   : regolith_streamtube_head
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : stream tube's head approximation  
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  The stream tube's head is always attached to the surface. 
!
!**********************************************************************************************************************************
subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots,age_collector,meltinejecta,totvol,distvol)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_head
   implicit none
   ! arguemnts
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(inout) :: meltinejecta,totvol
   real(DP),intent(in) :: deltar
   real(DP),intent(inout) :: totmare,tots
   real(SP),dimension(:),intent(inout) :: age_collector
   real(SP),dimension(:),intent(inout) :: distvol
 
   ! internal variables
   !type(regolisttype),pointer :: current
   type(regodatatype),dimension(:),allocatable :: current
   real(DP),parameter :: vratio = sqrt(2.0_DP)/2.0_DP ! Unfortunately, the approximate function that is used to get the size of a stream
                                                      ! tube with a constraint of CTEM's ejecta blanket thickness is slightly different
                                                      ! from the analytical function that we use here to approximate the stream tube's 
                                                      ! head intersected with underlying layers, about 30% of volume difference. 
   real(DP) :: z,zstart,zend,zmin,zmax
   real(DP) :: tothead,totmarehead,marehead,vhead,vsgly
   integer(I4B) :: N

   ! melt collector
   real(DP) :: recyratio
   real(DP) :: headmeltvol

   !current => surfi%regolayer
   !current = surfi%regolayer
   allocate(current,source=surfi%regolayer)
   N = size(current)
   z = current(N)%thickness
   vsgly = vratio * PI * deltar**3
   tothead = 0._DP
   totmarehead = 0._DP
   vhead = 0._DP
   marehead = 0.0_DP
   zstart = 0._DP
   zend = z
   zmin = zstart 
   zmax = 2.0 * deltar

   if (zend >= zmax) then ! Stream tube's head is inside the 1st layer.
      tots = tots + vsgly
      totmare = totmare + vsgly * current(N)%comp
      recyratio = vsgly / (user%pix**2) /current(N)%thickness
      age_collector(:) = age_collector(:) + current(N)%age(:) * recyratio
      headmeltvol = current(N)%meltfrac * recyratio * vsgly
      meltinejecta = meltinejecta + headmeltvol
      distvol(:) = distvol + (current(N)%meltdist(:)*recyratio*vsgly)
      totvol = totvol + vsgly
   else ! head is not intersected with layers. 

      do
         ! if (.not. associated(current%next)) exit
         
         if (zend < zmax) then 
            vhead = regolith_circle_sector_func(deltar,zstart,zend)
            tothead = tothead + vhead * vratio 
            totmarehead = totmarehead + vhead * vratio * current(N)%comp
            recyratio = vhead * vratio / (user%pix**2) / current(N)%thickness
            age_collector(:) = age_collector(:) + current(N)%age(:) * recyratio
            headmeltvol = current(N)%meltfrac * recyratio * tothead
            meltinejecta = meltinejecta + headmeltvol
            distvol(:) = distvol + (current(N)%meltdist(:)*recyratio*tothead)
            totvol = totvol + tothead
            !current => current%next
            N = N - 1
            z = z + current(N)%thickness
            zstart = zend
            zend = z
         else 
            totmarehead = totmarehead + (vsgly-tothead) * current(N)%comp
            tothead = vsgly
            recyratio = (vsgly - tothead) / (user%pix**2) / current(N)%thickness
            age_collector(:) = age_collector(:) + current(N)%age(:) * recyratio
            headmeltvol = current(N)%meltfrac * recyratio * tothead
            meltinejecta = meltinejecta + headmeltvol
            distvol(:) = distvol + (current(N)%meltdist(:)*recyratio*tothead)
            totvol = totvol + tothead
            exit
         end if
      end do


      tots = tots + tothead
      totmare = totmare + totmarehead

   end if

   return
end subroutine regolith_streamtube_head
