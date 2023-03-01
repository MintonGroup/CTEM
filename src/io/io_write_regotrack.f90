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
   integer(I4B), parameter :: FEJMF = 16
   integer(I4B), parameter :: FMF = 17
   integer(I4B), parameter :: FDF = 18
   !type(regolisttype),pointer :: current => null()
   type(regodatatype),dimension(:),allocatable :: current
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num
   real(DP),dimension(:),allocatable :: meltfrac, thickness, comp, ejm, ejmf, meltvolume
   real(SP),dimension(:,:),allocatable :: age, meltdist, distvol
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
   open(FMD,file=MDFILE,status='replace',form='unformatted')
   open(FEJM,file=EJMFILE,status='replace',form='unformatted')
   open(FMF,file=MELTFRACFILE,status='replace',form='unformatted')
   open(FDF,file=DISTFRACFILE,status='replace',form='unformatted')
   open(FEJMF,file=EJMFFILE,status='replace',form='unformatted')

   ! First pass to get stack numbers
   stacks_num(:,:) = 0
   do j=1,user%gridsize
      do i=1,user%gridsize
         !current => surf(i,j)%regolayer
         allocate(current,source=surf(i,j)%regolayer)
         stacks_num(i,j) = size(current)
         deallocate(current)
         ! do 
         !    if (.not. associated(current)) exit ! We've reached the bottom of the linked list
         !    stacks_num(i,j) = stacks_num(i,j) + 1
         !    current => current%next
         ! end do
      end do 
   end do

   ! Second pass to get data and save it
   do j=1,user%gridsize
      do i=1,user%gridsize
         !current => surf(i,j)%regolayer
         N = stacks_num(i,j)
         allocate(meltvolume(N),thickness(N),comp(N),age(MAXAGEBINS,N),distvol(domain%rcnum,N),ejm(N),ejmf(N),&
            meltfrac(N),meltdist(domain%rcnum,N))
         allocate(current,source=surf(i,j)%regolayer)
         do k=1,N
            meltfrac(k) = current(k)%meltfrac
            ! thickness(k) = current%regodata%thickness
            ! comp(k) = current%regodata%comp
            ! age(:,k) = current%regodata%age(:)
            ! current => current%next
            meltvolume(k) = current(k)%meltvolume
            thickness(k) = current(k)%thickness
            comp(k) = current(k)%comp
            age(:,k) = current(k)%age(:)
            !write(*,*) i, j
            distvol(:,k) = current(k)%distvol(:)
            meltdist(:,k) = current(k)%meltdist(:)
            ejm(k) = current(k)%ejm
            ejmf(k) = current(k)%ejmf

         end do
         deallocate(current)
         write(FMELT) meltvolume(:)
         write(FREGO) thickness(:)
         write(FCOMP) comp(:)
         write(FAGE) age(:,:)
         write(FMD) distvol(:,:)
         write(FEJM) ejm(:)
         write(FEJMF) ejmf(:)
         write(FMF) meltfrac(:)
         write(FDF) meltdist(:,:)
         deallocate(meltvolume,thickness,comp,age,distvol,ejm,ejmf,meltfrac,meltdist)
      end do 
   end do
   close(FMELT)
   close(FREGO)
   close(FCOMP)
   close(FAGE)
   close(FMD)
   close(FEJM)
   close(FEJMF)
   close(FMF)
   close(FDF)

   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUN,file=STACKNUMFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) stacks_num
   close(LUN)

   return
end subroutine io_write_regotrack