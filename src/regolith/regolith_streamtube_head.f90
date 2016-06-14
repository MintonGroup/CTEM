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
subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots)
!subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots,turnover)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_head
   implicit none
   ! arguemnts
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(in) :: deltar
   real(DP),intent(inout) :: totmare,tots

   ! internal variables
   type(regolisttype),pointer :: current
   real(DP),parameter :: vratio = sqrt(2.0_DP)/2.0_DP ! Unfortunately, the approximate function that is used to get the size of a stream
                                                      ! tube with a constraint of CTEM's ejecta blanket thickness is slightly different
                                                      ! from the analytical function that we use here to approximate the stream tube's 
                                                      ! head intersected with underlying layers, about 30% of volume difference. 
   real(DP) :: z,zstart,zend,zmin,zmax
   real(DP) :: tothead,totmarehead,marehead,vhead,vsgly

   ! * Mixing
   real(DP) :: zmix

   current => surfi%regolayer
   z = current%regodata%thickness
   zmix = z
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
      totmare = totmare + vsgly * current%regodata%comp
      !write(*,*) '0',zstart,zend,zmin,zmax,current%comp,totmare/2500.0,tots/2500.0
   else ! head is not intersected with layers. 

   do
      if (.not. associated(current%next)) exit
      
      if (zend < zmax) then 
         vhead = regolith_circle_sector_func(deltar,zstart,zend)
         tothead = tothead + vhead * vratio 
         totmarehead = totmarehead + vhead * vratio * current%regodata%comp
         !write(*,*) '1',zstart,zend,current%comp,vhead*vratio/2500.0,totmarehead/2500.0,tothead/2500.0
         current => current%next
         z = z + current%regodata%thickness
         zstart = zend
         zend = z
      else 
         totmarehead = totmarehead + (vsgly-tothead) * current%regodata%comp
         tothead = vsgly
         !write(*,*) '2',zstart,zend,current%comp,(vsgly-tothead)/2500.0,totmarehead/2500.0,tothead/2500.0
         exit
      end if
   end do

   tots = tots + tothead
   totmare = totmare + totmarehead

   end if

   return
end subroutine regolith_streamtube_head
