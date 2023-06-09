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
subroutine regolith_mix(user,surfi,mixing_depth,domain)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
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
   call util_traverse_pop_array(user,surfi%regolayer,mixing_depth,poppedarray)

   newlayer%thickness = 0.0_DP
   newlayer%comp      = 0.0_DP
   newlayer%age(:)    = 0.0_DP
   allocate(newlayer%distvol(1+domain%rcnum))
   newlayer%distvol(:) = 0.0_DP
   newlayer%ejm       = 0.0_DP
   newlayer%meltvolume = 0.0_DP
   newlayer%totvolume = 0.0_DP

   !poppedlist => poppedlist_top
   !do while(associated(poppedlist%next))
   N = size(poppedarray)
   do i = N,1,-1
      newlayer%thickness = newlayer%thickness + poppedarray(i)%thickness
      newlayer%comp      = newlayer%comp + poppedarray(i)%thickness * poppedarray(i)%comp       
      newlayer%age(:)    = newlayer%age(:) + poppedarray(i)%age(:)
      newlayer%distvol(:) = newlayer%distvol(:) + poppedarray(i)%distvol(:)
      newlayer%ejm       = newlayer%ejm + poppedarray(i)%ejm
      newlayer%meltvolume = newlayer%meltvolume + poppedarray(i)%meltvolume
   end do

   ! Get average values of composition and melt fraction
   newlayer%comp = newlayer%comp / newlayer%thickness 

   newlayer%totvolume = newlayer%thickness * user%pix * user%pix
   
   call util_push_array(surfi%regolayer, newlayer)
   !call util_destroy_list(poppedlist_top)


   ! do i = N,1,-1
   !    if (abs(surfi%regolayer(i)%meltvolume - sum(surfi%regolayer(i)%distvol) > 1e-5)) then
   !    write(*,*) "melt array =/= melt value!", domain%nqmc, domain%currentqmc, abs(surfi%regolayer(i)%meltvolume - sum(surfi%regolayer(i)%distvol))
   !    end if
   ! end do

   return
end subroutine regolith_mix
