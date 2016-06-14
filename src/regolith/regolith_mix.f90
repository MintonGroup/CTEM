!**********************************************************************************************************************************
!
!  Unit Name   : regolith_mix
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : mixing operation in the push-pop system        
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : surf : Surface expression matrix
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_mix(surfi,mixing_depth)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surfi
   real(DP), intent(in) :: mixing_depth

   ! Internal variables
   type(regolayertype) :: newlayer
   type(regolisttype),pointer :: poppedlist => null()

   !===============================================
   ! Add up all layers' info until a desired depth
   !===============================================          
   call util_traverse_pop(surfi%regolayer,mixing_depth,poppedlist)

   newlayer%thickness = 0.0_DP
   newlayer%comp = 0.0_DP
   newlayer%meltfrac = 0.0_DP

   do while(associated(poppedlist))
      newlayer%thickness = newlayer%thickness + poppedlist%regodata%thickness
      newlayer%comp = newlayer%comp + poppedlist%regodata%thickness * poppedlist%regodata%comp       
      newlayer%meltfrac = newlayer%meltfrac + poppedlist%regodata%thickness * poppedlist%regodata%meltfrac
      poppedlist => poppedlist%next
   end do

   ! Get average values of composition and melt fraction
   newlayer%comp = newlayer%comp / newlayer%thickness 
   newlayer%meltfac = newlayer%meltfrac / newlayer%thickness 
   
   call util_push(surf, newlayer)
   call util_destroy_list(poppedlist)

   return
end subroutine regolith_mix
