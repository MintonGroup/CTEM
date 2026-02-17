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
subroutine io_read_regotrack(user,surf,domain)
   use module_globals
   use module_util
   use module_io, EXCEPT_THIS_ONE => io_read_regotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(domaintype),intent(in)    :: domain

   ! Internals
   integer(I4B), parameter :: LUN=7
   integer(I4B), parameter :: FMELT = 10
   integer(I4B), parameter :: FREGO = 11
   integer(I4B), parameter :: FCOMP = 12
   integer(I4B), parameter :: FAGE = 13
   integer(I4B), parameter :: FMD = 14
   integer(I4B), parameter :: FEJM = 15
   ! real(DP),dimension(user%gridsize,user%gridsize) :: regotop,melt,comp,ejm,ejmf,meltfrac
   ! real(SP),dimension(user%gridsize,user%gridsize,domain%rcnum) :: meltdist, distfrac
   ! real(SP),dimension(user%gridsize,user%gridsize,MAXAGEBINS) :: age
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num 
   real(DP),dimension(:),allocatable :: regotop,melt,comp,ejm,thickness,meltvolume,agei
   real(SP),dimension(:,:),allocatable :: age, distvol
   !real(DP), dimension(:), allocatable :: regotopi,melti,compi,agei,dfi,ejmi,ejmfi,mdi,mfi
   type(regodatatype) :: newsurfi
   integer(I4B) :: ioerr,i,j,k,q,itmp,N
   integer(kind=8) :: recsize
   logical :: initstat
      
   ! Executable code

   ! Open a file for obtaining the number of stacks that is stored in each linked list
   ioerr = 0
   recsize = storage_size(itmp) * user%gridsize * user%gridsize / 8
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

   ! open(FCOMP,file=COMPFILE,status='old',form='unformatted',iostat=ioerr)
   ! if (ioerr/=0) then
   !     write(*,*) 'Error! Cannot read file ',trim(adjustl(COMPFILE))
   !     stop
   ! end if 

   open(FMELT,file=MELTFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(MELTFILE))
       stop
   end if   

   open(FEJM,file=EJMFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(EJMFILE))
       stop
   end if

   ! open(FMD,file=MDFILE,status='old',form='unformatted',iostat=ioerr)
   ! if (ioerr/=0) then 
   !     write(*,*) 'Error! Cannot read file ',trim(adjustl(MDFILE))
   !     stop
   ! end if  


   ! open(FAGE,file=AGEFILE,status='old',form='unformatted',iostat=ioerr)
   ! if (ioerr/=0) then
   !     write(*,*) 'Error! Cannot read file ',trim(adjustl(AGEFILE))
   !     stop
   ! end if

   ! Start pushing regolith thickness and melt fraction of each layer
   allocate(newsurfi%distvol(1+domain%rcnum))

   do j=1,user%gridsize
      do i=1,user%gridsize

         !call util_init_list(surf(i,j)%regolayer,initstat)
         !call util_init_array(user,surf(i,j)%regolayer,domain,initstat)
         N = stacks_num(i,j)
         allocate(meltvolume(N),thickness(N),comp(N),age(MAXAGEBINS,N),distvol(1+domain%rcnum,N),ejm(N))

         read(FMELT) meltvolume(:)
         read(FREGO) thickness(:)
         !read(FCOMP) comp(:)
         !read(FAGE) age(:,:)
         !read(FMD) distvol(:,:)
         read(FEJM) ejm(:)

         ! Temporary
         comp(:) = 0.0_DP
         age(:,:) = 0.0_DP
         distvol(:,:) = 0.0_DP
 
         allocate(agei(MAXAGEBINS * stacks_num(i,j)))

        
         do k=1,stacks_num(i,j)

            do q=1,MAXAGEBINS
               agei(MAXAGEBINS*k - (MAXAGEBINS-q)) = age(q,k)
            end do

         end do

         !do k=max(stacks_num(i,j),1),1,-1
         do k=1,max(stacks_num(i,j),1),1
            newsurfi%thickness = thickness(k)
            !newsurfi%comp = comp(k)
            newsurfi%meltvolume = meltvolume(k)
            newsurfi%ejm = ejm(k)
            newsurfi%totvolume = thickness(k) * user%gridsize * user%gridsize
            do q=1,MAXAGEBINS
               newsurfi%age(q) = agei(MAXAGEBINS*k-(MAXAGEBINS-q))
            end do
            do q=1,1+domain%rcnum
               newsurfi%distvol(q) = distvol(q,k)
            end do
            call util_push_array(surf(i,j)%regolayer,newsurfi)
         end do

         deallocate(meltvolume,thickness,comp,age,distvol,ejm,agei)

      end do
   end do
   close(FMELT)
   close(FREGO)
   ! close(FCOMP)
   ! close(FAGE)
   close(FEJM)
   !close(FMD)
   return
end subroutine io_read_regotrack
