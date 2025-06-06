!**********************************************************************************************************************************
!
!  Unit Name   : regolith_superdomain
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Softens the terrain under the ejecta  using a box filter model where the size of the box is proportional to the 
!                thickness of the ejecta  
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_superdomain(user,crater,domain,regolayer,ejdistribution,xpi,ypi,rm,depthb)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_superdomain
   implicit none

   ! Arguments
   type(usertype),intent(in)             :: user
   type(cratertype),intent(inout)        :: crater
   type(domaintype),intent(in)           :: domain
   type(regodatatype),dimension(:),allocatable,intent(inout)            :: regolayer
   real(DP),intent(in)                   :: ejdistribution
   integer(I4B),intent(in)               :: xpi, ypi
   real(DP),intent(in)                   :: rm 
   real(DP),intent(in)                   :: depthb

   ! Internal variables
   real(DP)                     :: xp, yp
   real(DP)                     :: lrad, vej
   real(DP)                     :: ebh
   real(DP)                     :: erad
   real(DP)                     :: cvpg
   real(DP)                     :: deltar, xmints, melt
   type(regodatatype)           :: newlayer   
   integer(I2B)                 :: n_age

   ! Interpret landing distance from a super domain crater to ray's deposits in
   ! local domain
   xp          = dble(xpi) * user%pix
   yp          = dble(ypi) * user%pix
   lrad        = sqrt((crater%xl - xp)**2 + (crater%yl - yp)**2)
   ebh         = 0.14_DP * crater%frad**(0.74_DP) * (lrad / crater%frad)**(-3.0_DP)   
   ebh         = ebh * ejdistribution
   ! Estimate radial position of a stream tube
   ! Assume streamlines' maximum depth is within regolith. 
   cvpg        = sqrt(2.0_DP) / 0.85_DP * ( user%mu_b / (user%mu_b + 1) )
   erad        = cvpg**(user%mu_b) * (lrad**(-0.5_DP * user%mu_b)) * (crater%rad)**(user%mu_b * 0.5_DP + 1.0_DP)
   vej         = cvpg * sqrt(user%gaccel * crater%grad) * (erad / crater%grad)**(-1.0_DP / user%mu_b) !equation 18 in Richardson 2009
   lrad        = ( vej **2 ) / user%gaccel !assume ejection angle is 45 degree.
   call regolith_melt_glass(user,crater,domain,ebh,rm,erad,lrad,deltar,newlayer,xmints,melt) 
   if (.not. allocated(newlayer%regotemp)) then
      allocate(newlayer%regotemp(1,1))
      newlayer%regotemp(1,1) = 0.0_SP
   end if
   if (.not. allocated(newlayer%regotime)) then
      allocate(newlayer%regotime(1,1))
      newlayer%regotime(1,1) = 0.0_SP
   end if
   call util_push_array(regolayer,newlayer)

   return
end subroutine regolith_superdomain
