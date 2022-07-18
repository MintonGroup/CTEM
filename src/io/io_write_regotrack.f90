!**********************************************************************************************************************************
!
!  Unit Name   : io_write_regotrack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for terrain grids for regolith tracking        
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
subroutine io_write_regotrack(user,surf)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_regotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf

   ! Regotrack Internals
   integer(I4B) :: i,j,k
   integer(I4B), parameter :: LUN = 7
   integer(I4B), parameter :: FMELT = 10
   integer(I4B), parameter :: FREGO = 11
   integer(I4B), parameter :: FCOMP = 12
   integer(I4B), parameter :: FAGE = 13
   type(regolisttype),pointer :: current => null()
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num
   real(DP),dimension(:),allocatable :: meltfrac, thickness, comp
   real(SP),dimension(:,:),allocatable :: age
   integer(kind=8) :: recsize
   real(DP) :: dtmp
   real(SP) :: stmp
   integer(I4B) :: itmp, N
   real(DP),dimension(user%gridsize,user%gridsize) :: comptop, rego
   real(DP),dimension(:),allocatable :: marehisto

   ! Executable code
   open(FMELT,file=MELTFILE,status='replace',form='unformatted')
   open(FREGO,file=REGOFILE,status='replace',form='unformatted')
   open(FCOMP,file=COMPFILE,status='replace',form='unformatted')
   open(FAGE,file=AGEFILE,status='replace',form='unformatted')

   ! First pass to get stack numbers
   stacks_num(:,:) = 0
   do j=1,user%gridsize
      do i=1,user%gridsize
         current => surf(i,j)%regolayer
         do 
            if (.not. associated(current)) exit ! We've reached the bottom of the linked list
            stacks_num(i,j) = stacks_num(i,j) + 1
            current => current%next
         end do
      end do 
   end do

   ! Second pass to get data and save it
   do j=1,user%gridsize
      do i=1,user%gridsize
         current => surf(i,j)%regolayer
         N = stacks_num(i,j)
         allocate(meltfrac(N),thickness(N),comp(N),age(MAXAGEBINS,N))
         do k=1,N
            meltfrac(k) = current%regodata%meltfrac
            thickness(k) = current%regodata%thickness
            comp(k) = current%regodata%comp
            age(:,k) = current%regodata%age(:)
            current => current%next
         end do
         write(FMELT) meltfrac(:)
         write(FREGO) thickness(:)
         write(FCOMP) comp(:)
         write(FAGE) age(:,:)
         deallocate(meltfrac,thickness,comp,age)
      end do 
   end do
   close(FMELT)
   close(FREGO)
   close(FCOMP)
   close(FAGE)

   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUN,file=STACKNUMFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) stacks_num
   close(LUN)

   return
end subroutine io_write_regotrack
