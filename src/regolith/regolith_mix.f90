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
subroutine regolith_mix(surf,d)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surf
   real(DP), intent(in) :: d

   ! Internal variables
   real(DP) :: z, zmare, ztot, z0
   type(regolayertype) :: newlayer

   z = surf%regolayer%thickness 
   zmare = 0._DP
   ztot  = 0._DP
   z0    = 0._DP

   !===============================================
   ! Add up all layers' info until a desired depth
   !===============================================          
   do while ( (associated(surf%regolayer%next)) .and. (z<d) ) 
            ztot  = ztot  + (z - z0)
            zmare = zmare + (z - z0) * surf%regolayer%comp
            call util_pop(surf)
            z0    = z
            z     = z + surf%regolayer%thickness
   end do

   ztot  = ztot  + (d - z0)
   zmare = zmare + (d - z0) * surf%regolayer%comp
   call regolith_traverse_pop(-1.0_DP * (d - z0), surf)
   newlayer%thickness = ztot 
   newlayer%comp      = zmare / ztot
   newlayer%meltfrac  = surf%regolayer%meltfrac
   call util_push(surf, newlayer)

end subroutine regolith_mix
