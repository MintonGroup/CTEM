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
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surfi
   real(DP), intent(in) :: mixing_depth

   ! Internal variables
   type(regolayertype) :: newlayer

   !===============================================
   ! Add up all layers' info until a desired depth
   !===============================================          
   call regolith_traverse_pop(mixing_depth, surfi,newlayer)
   call regolith_push(surf, newlayer)

end subroutine regolith_mix
