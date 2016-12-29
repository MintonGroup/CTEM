!**********************************************************************************************************************************
!
!  Unit Name   : io_read_porotrack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in files for pre-existing terrain grids for porosity files.
!
!  Input
!    Arguments : usr
!
!  Output
!    Arguments : surf : Surface expression matrix
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine io_read_porotrack(user,surf)
   use module_globals
   use module_util
   use module_io, EXCEPT_THIS_ONE => io_read_porotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(out) :: surf

   ! Internals
   integer(I4B), parameter :: LUN=7
   integer(I4B), parameter :: LUP=9
   integer(I4B), parameter :: LUC=10
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num 
   integer(I4B) :: ioerr, i, j, k, itmp, allocstat
   integer(kind=8) :: recsize
   type(regodatatype) :: newsurfi
   type(regolisttype), pointer :: current => null() 
   real(DP), dimension(:),allocatable :: depthk, porok
   logical :: initstat 
      
   ! Executable code

   ! Open a file for obtaining the number of stacks that is stored in each linked list
   ioerr = 0
   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUP,file=STACKPORFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(STACKPORFILE))
       stop
   end if
   read(LUP,rec=1) stacks_num
   close(LUP)

   ! Open files for depth. 
   ! Here, opening the unformatted file WITHOUT the record at the beginning and end of each value. 
   ! Follow the style in io_write_porotrack
   open(LUN,file=DEPTHFILE,status='old',form='unformatted',access='stream',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(DEPTHFILE))
       stop
   end if

   ! Open files for porosity. 
   ! Here, opening the unformatted file WITHOUT the record at the beginning and end of each value. 
   ! Follow the style in io_write_porotrack   
   open(LUC,file=POROFILE,status='old',form='unformatted',access='stream',iostat=ioerr)
   if (ioerr/=0) then
       write(*,*) 'Error! Cannot read file ',trim(adjustl(POROFILE))
       stop
   end if 

   do j=1,user%gridsize
      do i=1,user%gridsize

			! allocate depthk and porok, 
			! these are the layers at cell (i,j). 
         allocate(depthk(stacks_num(i,j)))
         allocate(porok(stacks_num(i,j)))
      
      	! from stacks_num to 1. 
         do k = stacks_num(i,j), 1, -1
            read(LUN) depthk(k)    
            read(LUC) porok(k)
         end do
         
         ! if surf(i,j)%porolayer is not allocated, allocate. 
         ! and add the deepest layer information on the top of the linked list layer. 
			if (.not. associated(surf(i,j)%porolayer)) then
				allocate(surf(i,j)%porolayer, STAT=allocstat)
				surf(i,j)%porolayer%regodata%depth    = depthk(1)
				surf(i,j)%porolayer%regodata%porosity = porok(1)			
			end if
			
			! if stacks_num is larger than 1, then keep pushing a new layer onto the old layer. 
         do k= 2, stacks_num(i,j)
            newsurfi%depth    = depthk(k)
            newsurfi%porosity = porok(k)
            call util_push(surf(i,j)%porolayer,newsurfi)
         end do

         deallocate(depthk, porok)

      end do
   end do
   close(LUN)
   close(LUC)
   return
end subroutine io_read_porotrack
