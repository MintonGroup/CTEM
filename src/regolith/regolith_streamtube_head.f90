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
subroutine regolith_streamtube_head(user,surfi,deltar,newlayer,eradi,rm)
!subroutine regolith_streamtube_head(user,surfi,deltar,totmare,tots,turnover)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_head
   implicit none
   ! arguemnts
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(in) :: deltar
   type(regodatatype),intent(inout) :: newlayer
   real(DP),intent(in) :: eradi,rm

   ! internal variables
   type(regolisttype),pointer :: current
   real(DP),parameter :: vratio = sqrt(2.0_DP)/2.0_DP ! Unfortunately, the approximate function that is used to get the size of a stream
                                                      ! tube with a constraint of CTEM's ejecta blanket thickness is slightly different
                                                      ! from the analytical function that we use here to approximate the stream tube's 
                                                      ! head intersected with underlying layers, about 30% of volume difference. 
   real(DP) :: z,zstart,zend,zmin,zmax
   real(DP) :: headtot,headcomp,headmelt,vhead,vsgly

   current => surfi%regolayer
   z = current%regodata%thickness
   vsgly = vratio * PI * deltar**3
   zstart = 0._DP
   zend = z
   zmin = zstart 
   zmax = 2.0 * deltar

   if (zend >= zmax) then ! Stream tube's head is inside the 1st layer.

      newlayer%thickness = newlayer%thickness + vsgly
      newlayer%comp      = newlayer%comp + vsgly * current%regodata%comp

   else ! head is not intersected with layers. 

   headtot = 0._DP
   headcomp = 0._DP
   headmelt = 0._DP
   vhead = 0._DP

   do
      if (.not. associated(current%next)) exit
      
      if (zend < zmax) then 
         vhead = regolith_circle_sector_func(deltar,zstart,zend)
         headtot = headtot + vhead * vratio 
         headcomp = headcomp + vhead * vratio * current%regodata%comp
         current => current%next
         z = z + current%regodata%thickness
         zstart = zend
         zend = z
      else 
         headcomp = headcomp + (vsgly-headtot) * current%regodata%comp
         headtot = vsgly
         !write(*,*) '2',zstart,zend,current%comp,(vsgly-tothead)/2500.0,totmarehead/2500.0,tothead/2500.0
         exit
      end if
   end do

   newlayer%thickness = newlayer%thickness + headtot
   newlayer%comp      = newlayer%comp      + headcomp

   end if

   return
end subroutine regolith_streamtube_head
