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
subroutine regolith_mix(surfi,mixing_depth,domain)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(surftype),intent(inout) :: surfi
   real(DP), intent(in) :: mixing_depth
   type(domaintype),intent(in) :: domain

   ! Internal variables
   type(regodatatype) :: newlayer
   !type(regolisttype),pointer :: poppedlist,poppedlist_top
   type(regodatatype),dimension(:),allocatable :: poppedarray
   integer(I4B) :: i, j, N

   !===============================================
   ! Add up all layers' info until a desired depth
   !=============================================== 
   ! !test code to create a situation for a breakpoint, since vscode debugger won't recognize the conditional breakpoint
   ! if(domain%currentqmc .eqv. .true.) then
   !    j = 0
   ! end if     
   call util_traverse_pop_array(surfi%regolayer,mixing_depth,poppedarray)

   newlayer%thickness = 0.0_DP
   newlayer%comp      = 0.0_DP
   newlayer%meltfrac  = 0.0_DP
   newlayer%age(:)    = 0.0_DP
   allocate(newlayer%meltdist(domain%rcnum))
   newlayer%meltdist(:) = 0.0_SP

   !poppedlist => poppedlist_top
   !do while(associated(poppedlist%next))
   N = size(poppedarray)
   do i = N,1,-1
      newlayer%thickness = newlayer%thickness + poppedarray(i)%thickness
      newlayer%comp      = newlayer%comp + poppedarray(i)%thickness * poppedarray(i)%comp       
      newlayer%meltfrac  = newlayer%meltfrac + poppedarray(i)%thickness * poppedarray(i)%meltfrac
      newlayer%age(:)    = newlayer%age(:) + poppedarray(i)%age(:)
      newlayer%meltdist(:) = newlayer%meltdist(:) + poppedarray(i)%thickness * poppedarray(i)%meltdist(:)
      ! do j = 1,domain%rcnum !testing a loop here since the array operation resulted in a segfault
      !    newlayer%meltdist(j) = newlayer%meltdist(j) + poppedarray(i)%thickness * poppedarray(i)%meltdist(j)
      ! end do
      ! !poppedlist => poppedlist%next
   end do

   ! Get average values of composition and melt fraction
   newlayer%comp = newlayer%comp / newlayer%thickness 
   newlayer%meltfrac = newlayer%meltfrac / newlayer%thickness 
   newlayer%meltdist(:) = newlayer%meltdist(:) / newlayer%thickness
   
   call util_push_array(surfi%regolayer, newlayer)
   !call util_destroy_list(poppedlist_top)

   return
end subroutine regolith_mix
