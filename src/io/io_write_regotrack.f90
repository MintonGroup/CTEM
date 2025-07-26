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
subroutine io_write_regotrack(user,surf,domain)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_regotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   type(domaintype),intent(in) :: domain

   ! Regotrack Internals
   integer(I4B) :: i,j,k
   integer(I4B), parameter :: LUN = 7
   integer(I4B), parameter :: FMELT = 10
   integer(I4B), parameter :: FREGO = 11
   integer(I4B), parameter :: FCOMP = 12
   integer(I4B), parameter :: FAGE = 13
   integer(I4B), parameter :: FMD = 14
   integer(I4B), parameter :: FEJM = 15
   integer(I4B), parameter :: FRT = 16
   integer(I4B), parameter :: FT = 17
   integer(I4B), parameter :: FF = 18
   !type(regolisttype),pointer :: current => null()
   type(regodatatype),dimension(:),allocatable :: current
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num
   real(DP),dimension(:),allocatable :: thickness, comp, ejm, meltvolume
   real(SP),dimension(:,:),allocatable :: regotemp, regotime
   real(SP),dimension(:,:),allocatable :: age, distvol
   real(SP),dimension(:),allocatable :: frac
   integer(kind=8) :: recsize
   real(DP) :: dtmp
   real(SP) :: stmp
   integer(I4B) :: itmp, N, NT, t, NC, c, l, m
   real(DP),dimension(user%gridsize,user%gridsize) :: comptop, rego
   real(DP),dimension(:),allocatable :: marehisto

   ! Executable code
   open(FMELT,file=MELTFILE,status='replace',form='unformatted')
   open(FREGO,file=REGOFILE,status='replace',form='unformatted')
   open(FCOMP,file=COMPFILE,status='replace',form='unformatted')
   open(FAGE,file=AGEFILE,status='replace',form='unformatted')
   open(FMD,file=MDFILE,status='replace',form='unformatted')
   open(FEJM,file=EJMFILE,status='replace',form='unformatted')
   if (user%dothermal) then
      open(FRT,file=REGOTEMPFILE,status='replace',form='unformatted')
      open(FT,file=REGOTIMEFILE,status='replace',form='unformatted')
      open(FF,file=FRACFILE,status='replace',form='unformatted')
   end if

   ! First pass to get stack numbers
   stacks_num(:,:) = 0
   do j=1,user%gridsize
      do i=1,user%gridsize
         !current => surf(i,j)%regolayer
         allocate(current,source=surf(i,j)%regolayer)
         stacks_num(i,j) = size(current)
         deallocate(current)
      end do 
   end do

   ! Second pass to get data and save it
   do j=1,user%gridsize
      do i=1,user%gridsize
         !current => surf(i,j)%regolayer
         N = stacks_num(i,j)
         allocate(meltvolume(N),thickness(N),comp(N),age(MAXAGEBINS,N),distvol(1+domain%rcnum,N),ejm(N))
         allocate(current,source=surf(i,j)%regolayer)
         do k=1,N
            meltvolume(k) = current(k)%meltvolume
            thickness(k) = current(k)%thickness
            comp(k) = current(k)%comp
            age(:,k) = current(k)%age(:)
            !write(*,*) i, j
            distvol(:,k) = current(k)%distvol(:)
            ejm(k) = current(k)%ejm
            if (user%dothermal) then
               ! NT = size(current(k)%regotemp, 1)
               ! NC = size(current(k)%regotemp, 2)
               ! if(.not. allocated(regotemp)) allocate(regotemp(NC,NT))
               ! if(.not. allocated(regotime)) allocate(regotime(NC,NT))
               ! do t = 1, NT
               !    do c = 1, NC
               ! print *, 'regotemp bounds:', lbound(current(k)%regotemp), ubound(current(k)%regotemp)
               ! print *, 'regotemp size:', size(current(k)%regotemp)
               t = size(current(k)%regotemp,1)
               c = size(current(k)%regotemp,2)
               write(FRT) ((current(k)%regotemp(l,m), m=1,c), l=1,t)
               write(FT) ((current(k)%regotime(l,m), m=1,c), l=1,t)
               write(FF) current(k)%frac
               !    end do
               ! end do
            end if
         end do
         deallocate(current)
         write(FMELT) meltvolume(:)
         write(FREGO) thickness(:)
         write(FCOMP) comp(:)
         write(FAGE) age(:,:)
         write(FMD) distvol(:,:)
         write(FEJM) ejm(:)
         deallocate(meltvolume,thickness,comp,age,distvol,ejm)
         !if (user%dothermal) deallocate(regotemp,regotime)
      end do 
   end do
   close(FMELT)
   close(FREGO)
   close(FCOMP)
   close(FAGE)
   close(FMD)
   close(FEJM)
   if (user%dothermal) then
      close(FRT)
      close(FT)
      close(FF)
   end if

   recsize = storage_size(itmp) * user%gridsize * user%gridsize / 8
   open(LUN,file=STACKNUMFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) stacks_num
   close(LUN)

   return
end subroutine io_write_regotrack