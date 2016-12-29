!**********************************************************************************************************************************
!
!  Unit Name   : io_write_porotrack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for terrain grids for porosity tracking       
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
subroutine io_write_porotrack(user,surf)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_porotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf

   ! Regotrack Internals
   type(regolisttype), pointer :: current => null()
   integer(I4B) :: i,j,k
   integer(I4B) :: itmp
   integer(I4B), parameter :: LUM=8
   integer(I4B), parameter :: LUN=7
   integer(I4B), dimension(user%gridsize, user%gridsize) :: stacks_num
   integer(kind=8) :: recsize
   real(DP) :: dtmp
   real(DP) :: depth_ij, porosity_ij
   
   ! Executable code
   
   ! Here, we remove the record of "unformatted" values at the beginning and end of them. 
   ! To include the record, remove access = 'stream'
   open(LUM,file=DEPTHFILE,status='replace', form='unformatted', access = 'stream')
   open(LUN,file=POROFILE,status='replace', form='unformatted', access = 'stream')
   
   do j=1,user%gridsize
      do i=1,user%gridsize
         stacks_num(i,j) = 0
         current => surf(i,j)%porolayer
         ! output the porosity and the depth
         do 
         	if (.not. associated(current)) exit
         	stacks_num(i,j) = stacks_num(i,j) + 1
         	depth_ij    = current%regodata%depth
         	porosity_ij = current%regodata%porosity
         	
         	write(LUM) depth_ij
         	write(LUN) porosity_ij
         	current => current%next
         end do
      end do 
   end do
   close(LUN)
   close(LUM)

   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUN,file=STACKPORFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) stacks_num
   close(LUN)

   return
end subroutine io_write_porotrack
