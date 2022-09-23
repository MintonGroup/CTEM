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
   use module_util
   use module_io, EXCEPT_THIS_ONE => io_read_regotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(out) :: surf

   ! Internals
   integer(I4B), parameter :: LUN=7
   integer(I4B), parameter :: FMELT = 10
   integer(I4B), parameter :: FREGO = 11
   integer(I4B), parameter :: FCOMP = 12
   integer(I4B), parameter :: FAGE = 13
   real(DP),dimension(user%gridsize,user%gridsize) :: regotop,melt,comp
   real(SP),dimension(user%gridsize,user%gridsize,MAXAGEBINS) :: age
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num 
   real(DP), dimension(:), allocatable :: regotopi,melti,compi,agei
   type(regodatatype) :: newsurfi
   integer(I4B) :: ioerr,i,j,k,q,itmp
   integer(kind=8) :: recsize
   logical :: initstat 
      
   ! Executable code

   ! Open a file for obtaining the number of stacks that is stored in each linked list
   ioerr = 0
   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUN,file=STACKNUMFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(STACKNUMFILE))
       stop
   end if
   read(LUN,rec=1) stacks_num
   close(LUN)

   ! Open files for regolith thickness/melt fraction 
   open(FREGO,file=REGOFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(REGOFILE))
       stop
   end if

   open(FCOMP,file=COMPFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then
       write(*,*) 'Error! Cannot read file ',trim(adjustl(COMPFILE))
       stop
   end if 

   open(FMELT,file=MELTFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(MELTFILE))
       stop
   end if   

   open(FAGE,file=AGEFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then
       write(*,*) 'Error! Cannot read file ',trim(adjustl(MELTFILE))
       stop
   end if

   ! Start pushing regolith thickness and melt fraction of each layer in FILO manner

   do j=1,user%gridsize
      do i=1,user%gridsize

         call util_init_list(surf(i,j)%regolayer,initstat)

         allocate(regotopi(stacks_num(i,j)))
         allocate(compi(stacks_num(i,j)))
         allocate(melti(stacks_num(i,j)))    
         allocate(agei(MAXAGEBINS * stacks_num(i,j)))
        
         do k=1,stacks_num(i,j)
            read(FREGO) regotop(i,j)
            regotopi(k) = regotop(i,j)       
            read(FCOMP) comp(i,j)
            compi(k) = comp(i,j)
            read(FMELT) melt(i,j) 
            melti(k) = melt(i,j)
            read(FAGE) age(i,j,:)
            do q=1,MAXAGEBINS
               agei(MAXAGEBINS*k - (MAXAGEBINS-q)) = age(i,j,q)
            end do
         end do

         do k=max(stacks_num(i,j)-1,1),1,-1
            newsurfi%thickness = regotopi(k)
            newsurfi%comp = compi(k)
            newsurfi%meltfrac  = melti(k)
            do q=1,MAXAGEBINS
               newsurfi%age(q) = agei(MAXAGEBINS*k-(MAXAGEBINS-q))
            end do
            call util_push(surf(i,j)%regolayer,newsurfi)
         end do 

         deallocate(regotopi,compi,melti,agei)

      end do
   end do
   close(FMELT)
   close(FREGO)
   close(FCOMP)
   close(FAGE)
   return
end subroutine io_read_regotrack
