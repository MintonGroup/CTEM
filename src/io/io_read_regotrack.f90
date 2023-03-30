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
   type(surftype),dimension(:,:),intent(out) :: surf
   type(domaintype),intent(in)    :: domain

   ! Internals
   integer(I4B), parameter :: LUN=7
   integer(I4B), parameter :: FMELT = 10
   integer(I4B), parameter :: FREGO = 11
   integer(I4B), parameter :: FCOMP = 12
   integer(I4B), parameter :: FAGE = 13
   integer(I4B), parameter :: FMD = 14
   integer(I4B), parameter :: FEJM = 15
   integer(I4B), parameter :: FEJMF = 16
   integer(I4B), parameter :: FMF = 17
   integer(I4B), parameter :: FDF = 18
   ! real(DP),dimension(user%gridsize,user%gridsize) :: regotop,melt,comp,ejm,ejmf,meltfrac
   ! real(SP),dimension(user%gridsize,user%gridsize,domain%rcnum) :: meltdist, distfrac
   ! real(SP),dimension(user%gridsize,user%gridsize,MAXAGEBINS) :: age
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num 
   real(DP),dimension(:),allocatable :: regotop,melt,comp,ejm,ejmf,meltfrac,thickness,meltvolume
   real(SP),dimension(:,:),allocatable :: age, meltdist, distfrac, distvol
   real(DP), dimension(:), allocatable :: regotopi,melti,compi,agei,dfi,ejmi,ejmfi,mdi,mfi
   type(regodatatype) :: newsurfi
   integer(I4B) :: ioerr,i,j,k,q,itmp,N
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

   open(FDF,file=DISTFRACFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(DISTFRACFILE))
       stop
   end if  

   open(FEJM,file=EJMFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(EJMFILE))
       stop
   end if
   
   open(FEJMF,file=EJMFFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(EJMFFILE))
       stop
   end if  

   open(FMD,file=MDFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(MDFILE))
       stop
   end if  

   open(FMF,file=MELTFRACFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then 
       write(*,*) 'Error! Cannot read file ',trim(adjustl(MELTFRACFILE))
       stop
   end if  

   open(FAGE,file=AGEFILE,status='old',form='unformatted',iostat=ioerr)
   if (ioerr/=0) then
       write(*,*) 'Error! Cannot read file ',trim(adjustl(AGEFILE))
       stop
   end if

   ! Start pushing regolith thickness and melt fraction of each layer in FILO manner

   allocate(newsurfi%meltdist(domain%rcnum))
   allocate(newsurfi%distvol(domain%rcnum))

   do j=1,user%gridsize
      do i=1,user%gridsize

         !call util_init_list(surf(i,j)%regolayer,initstat)
         call util_init_array(user,surf(i,j)%regolayer,domain,initstat)
         N = stacks_num(i,j)
         allocate(meltvolume(N),thickness(N),comp(N),age(MAXAGEBINS,N),distvol(domain%rcnum,N),ejm(N),ejmf(N),&
            meltfrac(N),meltdist(domain%rcnum,N))

         read(FMELT) meltvolume(:)
         read(FREGO) thickness(:)
         read(FCOMP) comp(:)
         read(FAGE) age(:,:)
         read(FMD) distvol(:,:)
         read(FEJM) ejm(:)
         read(FEJMF) ejmf(:)
         read(FMF) meltfrac(:)
         read(FDF) meltdist(:,:)

         allocate(regotopi(stacks_num(i,j)))
         allocate(compi(stacks_num(i,j)))
         allocate(melti(stacks_num(i,j)))    
         allocate(agei(MAXAGEBINS * stacks_num(i,j)))
         allocate(ejmi(stacks_num(i,j)))
         allocate(ejmfi(stacks_num(i,j)))
         allocate(mfi(stacks_num(i,j)))
         allocate(dfi(domain%rcnum * stacks_num(i,j)))
         allocate(mdi(domain%rcnum * stacks_num(i,j)))

         ! allocate(melti(N),regotopi(N),compi(N),agei(MAXAGEBINS,N),mdi(domain%rcnum,N),ejmi(N),ejmfi(N),&
         !    mfi(N),mdi(domain%rcnum,N))
        
         do k=1,stacks_num(i,j)
         !    read(FREGO) regotop(i,j)
         !    regotopi(k) = regotop(i,j)       
         !    read(FCOMP) comp(i,j)
         !    compi(k) = comp(i,j)
         !    read(FMELT) melt(i,j) 
         !    melti(k) = melt(i,j)
         !    read(FMF) meltfrac(i,j)
         !    mfi(k) = meltfrac(i,j)
         !    read(FEJM) ejm(i,j)
         !    ejmi(k) = ejm(i,j)
         !    read(FEJMF) ejmf(i,j)
         !    ejmfi(k) = ejmf(i,j)
         !    read(FDF) distfrac(i,j,:)
         !    read(FMD) meltdist(i,j,:)
         !    do q=1,domain%rcnum
         !       dfi(q*k) = distfrac(i,j,q)
         !       mdi(q*k) = meltdist(i,j,q)
         !    end do
         !    read(FAGE) age(i,j,:)
         !    do q=1,MAXAGEBINS
         !       agei(MAXAGEBINS*k - (MAXAGEBINS-q)) = age(i,j,q)
         !    end do

            regotopi(k) = regotop(k)
            compi(k) = comp(k)
            ejmi(k) = ejm(k)
            ejmfi(k) = ejmf(k)
            mfi(k) = meltfrac(k)
            do q=1,domain%rcnum
               dfi(q*k) = distfrac(q,k) !or is it (k,q) ?
               mdi(q*k) = meltdist(q,k)
            end do
            do q=1,MAXAGEBINS
               agei(MAXAGEBINS*k - (MAXAGEBINS-q)) = age(q,k) !again, (k,q)?
            end do

         end do
            



         do k=max(stacks_num(i,j)-1,1),1,-1
            newsurfi%thickness = regotopi(k)
            newsurfi%comp = compi(k)
            newsurfi%meltfrac  = mfi(k)
            newsurfi%meltvolume = melti(k)
            newsurfi%ejm = ejmi(k)
            newsurfi%ejmf = ejmfi(k)
            newsurfi%totvolume = regotopi(k) * user%gridsize * user%gridsize
            do q=1,MAXAGEBINS
               newsurfi%age(q) = agei(MAXAGEBINS*k-(MAXAGEBINS-q))
            end do
            do q=1,domain%rcnum
               newsurfi%meltdist(q) = dfi(q*k)
               newsurfi%distvol(q) = mdi(q*k) !these two could be wrong based on the way the file is read; need to check
            end do
            call util_push_array(surf(i,j)%regolayer,newsurfi)
         end do 

         deallocate(regotopi,compi,melti,agei,dfi,ejmi,ejmfi,mdi,mfi)
         deallocate(meltvolume,thickness,comp,age,distvol,ejm,ejmf,meltfrac,meltdist)

      end do
   end do
   close(FMELT)
   close(FREGO)
   close(FCOMP)
   close(FAGE)
   close(FDF)
   close(FEJM)
   close(FEJMF)
   close(FMD)
   close(FMF)
   return
end subroutine io_read_regotrack
