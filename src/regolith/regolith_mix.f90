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
   integer(I4B) :: i, j, k, N, NT, X

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
   newlayer%age(:)    = 0.0_SP
   allocate(newlayer%distvol(1+domain%rcnum))
   newlayer%distvol(:) = 0.0_SP
   newlayer%ejm       = 0.0_DP
   newlayer%meltvolume = 0.0_DP
   newlayer%totvolume = 0.0_DP

   !poppedlist => poppedlist_top
   !do while(associated(poppedlist%next))
   N = size(poppedarray)
   NT = 0
   do i = N,1,-1
      newlayer%thickness = newlayer%thickness + poppedarray(i)%thickness
      newlayer%comp      = newlayer%comp + poppedarray(i)%thickness * poppedarray(i)%comp       
      newlayer%age(:)    = newlayer%age(:) + poppedarray(i)%age(:)
      newlayer%distvol(:) = newlayer%distvol(:) + poppedarray(i)%distvol(:)
      newlayer%ejm       = newlayer%ejm + poppedarray(i)%ejm
      newlayer%meltvolume = newlayer%meltvolume + poppedarray(i)%meltvolume
      X = size(poppedarray(i)%regotemp)
      if (X > NT) then 
         X = NT
      end if
   end do

   ! Get average values of composition and melt fraction
   newlayer%comp = newlayer%comp / newlayer%thickness 

   newlayer%totvolume = newlayer%thickness * user%pix * user%pix

   if (.not. allocated(newlayer%regotemp)) then
      allocate(newlayer%regotemp(1,1))
      newlayer%regotemp(1,1) = 0.0_SP
   end if
   if (.not. allocated(newlayer%regotime)) then
      allocate(newlayer%regotime(1,1))
      newlayer%regotime(1,1) = 0.0_SP
   end if
   if (.not. allocated(newlayer%frac)) then
      allocate(newlayer%frac(1))
      newlayer%frac = 0.0_SP
   end if
   
   do i=1,N
      X = size(poppedarray(i)%regotemp)
      do j=1,NT
         if (j > X) then
            newlayer%regotemp(i,j) = -1.0_SP !NoData value
         else
            newlayer%regotemp(i,j) = poppedarray(i)%regotemp(i,j)
         end if
      end do
   end do
   
   call util_push_array(surfi%regolayer, newlayer)


   return
end subroutine regolith_mix
