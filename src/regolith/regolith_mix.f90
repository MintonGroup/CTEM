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
   type(regodatatype) :: newlayer
   !type(regolisttype),pointer :: poppedlist,poppedlist_top
   type(regodatatype),dimension(:),allocatable :: poppedarray, poppedarray_top
   integer(I4B) :: N

   !===============================================
   ! Add up all layers' info until a desired depth
   !===============================================          
   call util_traverse_pop_array(surfi%regolayer,mixing_depth,poppedarray_top)

   newlayer%thickness = 0.0_DP
   newlayer%comp      = 0.0_DP
   newlayer%meltfrac  = 0.0_DP
   newlayer%age(:)    = 0.0_DP

   !poppedlist => poppedlist_top
   !do while(associated(poppedlist%next))
   N = size(poppedarray)
   do
      newlayer%thickness = newlayer%thickness + poppedarray(N)%thickness
      newlayer%comp      = newlayer%comp + poppedarray(N)%thickness * poppedarray(N)%comp       
      newlayer%meltfrac  = newlayer%meltfrac + poppedarray(N)%thickness * poppedarray(N)%meltfrac
      newlayer%age(:)    = newlayer%age(:) + poppedarray(N)%age(:)
      !poppedlist => poppedlist%next
      N = N - 1
   end do

   ! Get average values of composition and melt fraction
   newlayer%comp = newlayer%comp / newlayer%thickness 
   newlayer%meltfrac = newlayer%meltfrac / newlayer%thickness 
   
   call util_push_array(surfi%regolayer, newlayer)
   !call util_destroy_list(poppedlist_top)
   deallocate(poppedarray_top)

   return
end subroutine regolith_mix
