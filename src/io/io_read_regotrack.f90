!**********************************************************************************************************************************
!
!  Unit Name   : io_read_regotrack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in files for pre-existing terrain grids for regolith tracking material
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
subroutine io_read_regotrack(user,surf)
   use module_globals
   use module_regolith
   use module_io, EXCEPT_THIS_ONE => io_read_regotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(out) :: surf

   ! Internals
   integer(I4B), parameter :: LUN=7
   integer(I4B), parameter :: LUM=8
   integer(I4B), parameter :: LUP=9
   integer(I4B), parameter :: LUC=10
   real(DP),dimension(user%gridsize,user%gridsize) :: regotop,melt,comp!, mixdep, excavdep
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num!, mix_num 
   real(DP), dimension(:),allocatable :: regotopi,melti,compi
   type(regolayertype) :: newlayer
   integer(I4B) :: ioerr,i,j,k,itmp
   integer(kind=8) :: recsize
   real(DP) :: dtmp
   ! Checking if it reads a pre-existing terrain right
   !integer(kind=8) :: recsize
   !real(DP) :: dtmp
   !real(DP),dimension(:,:),allocatable :: surface
      
   ! Executable code

   ! Open a file for obtaining the number of stacks that is stored in each linked list
   ioerr = 0
   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUP,file=STACKNUMFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(STACKNUMFILE))
       stop
   end if
   read(LUP,rec=1) stacks_num
   close(LUP)

   ! Open files for regolith thickness/melt fraction 
   open(LUN,file=REGOFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(REGOFILE))
       stop
   end if

   open(LUC,file=COMPFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then
       write(*,*) 'Error! Cannot read file ',trim(adjustl(COMPFILE))
       stop
   end if 

   !open(LUM,file=MELTFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(MELTFILE))
       stop
   end if   

   ! Start pushing regolith thickness and melt fraction of each layer in FILO manner
   do j=1,user%gridsize
      do i=1,user%gridsize

         allocate(regotopi(stacks_num(i,j)))
         allocate(compi(stacks_num(i,j)))
         allocate(melti(stacks_num(i,j)))    
      
         do k=1,stacks_num(i,j)
            read(LUN) regotop(i,j)
            regotopi(k) = regotop(i,j)       
            read(LUC) comp(i,j)
            compi(k) = comp(i,j)
            read(LUM) melt(i,j) 
            melti(k) = melt(i,j)
            !write(*,*) i,j,k,regotopi(k)
         end do

         allocate(surf(i,j)%regolayer)
         nullify(surf(i,j)%regolayer%next)

         do k=stacks_num(i,j),1,-1
            newlayer%thickness = regotopi(k)
            newlayer%comp = compi(k)
            newlayer%meltfrac  = melti(k)
            call regolith_push(surf(i,j),newlayer)
            !write(*,*) i,j,k,surf(i,j)%regolayer%thickness
         end do 

         deallocate(regotopi,compi,melti)

      end do
   end do
   close(LUN)
   close(LUC)
   close(LUM)
   return
end subroutine io_read_regotrack
