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
   integer(I4B) :: i, j, k, N, NT, X, total_rows, row_start, rows, cols, r, c, N2
   real(DP) :: tot


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
   N2 = size(surfi%regolayer)
   NT = 0
   do i = N,1,-1
      newlayer%thickness = newlayer%thickness + poppedarray(i)%thickness
      newlayer%comp      = newlayer%comp + poppedarray(i)%thickness * poppedarray(i)%comp       
      newlayer%age(:)    = newlayer%age(:) + poppedarray(i)%age(:)
      newlayer%distvol(:) = newlayer%distvol(:) + poppedarray(i)%distvol(:)
      newlayer%ejm       = newlayer%ejm + poppedarray(i)%ejm
      newlayer%meltvolume = newlayer%meltvolume + poppedarray(i)%meltvolume
      if (user%dothermal) then
         X = size(poppedarray(i)%regotemp)
         if (X > NT) then 
            NT = X
         end if
      end if
   end do

   !if (user%dothermal) call regolith_combine_temperatures(newlayer,poppedarray,N)

   ! allocate(newlayer%regotemp(NT))
   ! allocate(newlayer%regotime(NT))
   ! newlayer%regotemp(:) = -1.0_DP
   ! newlayer%regotime(:) = -1.0_DP


   total_rows = 0
   do i = 1, N
       if (allocated(poppedarray(i)%regotemp)) then
           total_rows = total_rows + size(poppedarray(i)%regotemp, 1)
       end if
   end do

  !  deallocate(poppedarray(1)%regotemp,poppedarray(1)%regotime)
  !  allocate(poppedarray(1)%regotemp(total_rows, NT))
  !  allocate(poppedarray(1)%regotime(total_rows, NT))
   allocate(newlayer%regotemp(total_rows, NT))
   allocate(newlayer%regotime(total_rows, NT))
   allocate(newlayer%frac(total_rows))
   newlayer%regotemp = 0.0_SP  ! Fill with NoData initially
   newlayer%regotime = 0.0_SP

   row_start = 1
   do i = 1, N
       if (allocated(poppedarray(i)%regotemp)) then
           rows = size(poppedarray(i)%regotemp, 1)
           cols = size(poppedarray(i)%regotemp, 2)
           do r = 1, rows
              newlayer%frac(row_start + r - 1) = poppedarray(i)%thickness / mixing_depth
               do c = 1, cols
                   newlayer%regotemp(row_start + r - 1, c) = poppedarray(i)%regotemp(r, c)
                   newlayer%regotime(row_start + r - 1, c) = poppedarray(i)%regotime(r, c)
               end do
           end do
           row_start = row_start + rows
       end if
   end do

   tot = sum(newlayer%frac)
   newlayer%frac(:) = newlayer%frac(:) / tot



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
      newlayer%frac = 1.0_SP
  end if
   
   call util_push_array(surfi%regolayer, newlayer)
   !call util_destroy_list(poppedlist_top)

   !deallocate(regotemps,regotimes)


   ! do i = N,1,-1
   !    if (abs(surfi%regolayer(i)%meltvolume - sum(surfi%regolayer(i)%distvol) > 1e-5)) then
   !    write(*,*) "melt array =/= melt value!", domain%nqmc, domain%currentqmc, abs(surfi%regolayer(i)%meltvolume - sum(surfi%regolayer(i)%distvol))
   !    end if
   ! end do

   return
end subroutine regolith_mix
