!**********************************************************************************************************************************
!
!  Unit Name   : crust_cta_generate
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Generates the crustal thin area below a big impact caused by  mantle uplift
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
subroutine crust_cta_generate(user,surfi,crater,domain,lradsq,depth,mdepth,trans_depth)
   use module_globals
   use module_util
   use module_crust, EXCEPT_THIS_ONE => crust_cta_generate
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: lradsq
   real(DP),intent(in) :: depth
   real(DP),intent(in) :: mdepth
   real(DP),intent(in) :: trans_depth

   ! Internal variables
   real(DP),parameter :: mindepth = 5000.0_DP ! Minimum crustal thinning (set arbitrarily to 5 km at the moment)
   real(DP) :: newdepth,parab

   ! Executable code

   parab = crater%parab * crater%rad/crater%frad
   newdepth = max(mindepth,mdepth - (trans_depth - parab*lradsq))
   surfi%mantle = surfi%dem - newdepth

   return
end subroutine crust_cta_generate

