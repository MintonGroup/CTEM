!**********************************************************************************************************************************
!
!  Unit Name   : io_read_surf
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in files for pre-existing terrain grids
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
subroutine io_read_surf(user,surf,domain)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_read_surf
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(out) :: surf
   type(domaintype),intent(in)    :: domain

   ! Internals
   integer(I4B) :: ioerr,i
   integer(kind=8) :: recsize
   integer(I4B), parameter :: LUN=7
   real(DP) :: dtmp
   real(SP) :: stmp
   integer(I2B) :: itmp

   ! Executable code

   ! Read in matrix files
   recsize = sizeof(dtmp) * user%gridsize * user%gridsize
   !write(*,*) 'recsize = ',recsize
   !write(*,*) 'alternate = ', sizeof(dtmp) * user%gridsize * user%gridsize
   !read(*,*) 
   open(LUN,file=DEMFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then
      write(*,*) 'Error! Cannot read file ',trim(adjustl(DEMFILE))
      stop
   end if
   read(LUN,rec=1) surf%dem
   close(LUN)


   open(LUN,file=EJCOVFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then
      write(*,*) 'Error! Cannot read file ',trim(adjustl(EJCOVFILE))
      stop
   end if
   read(LUN,rec=1) surf%ejcov
   close(LUN)

   open(LUN,file=DIAMFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then
      write(*,*) 'Error! Cannot read file ',trim(adjustl(DIAMFILE))
      stop
   end if
   do i=1,user%numlayers 
      read(LUN,rec=i) surf%diam(i)
   end do
   close(LUN)

   open(LUN,file=TIMEFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then
      write(*,*) 'Error! Cannot read file ',trim(adjustl(TIMEFILE))
      stop
   end if
   do i=1,user%numlayers 
      read(LUN,rec=i) surf%timestamp(i)
   end do
   close(LUN)

   recsize = sizeof(stmp) * user%gridsize * user%gridsize
   open(LUN,file=POSFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   if (ioerr/=0) then
      write(*,*) 'Error! Cannot read file ',trim(adjustl(POSFILE))
      stop
   end if
   do i=1,user%numlayers 
      read(LUN,rec=2*i-1) surf%xl(i)
      read(LUN,rec=2*i) surf%yl(i)
   end do
   close(LUN)

   if (user%doregotrack) call io_read_regotrack(user,surf,domain)
   
   !if (user%docrustal_thinning) then
   !   recsize=sizeof(itmp)*user%gridsize*user%gridsize
   !   open(LUN,file=THICKFILE,status='old',form='unformatted',recl=recsize,access='direct',iostat=ioerr)
   !   if (ioerr/=0) then
   !      write(*,*) 'Error! Cannot read file ',trim(adjustl(THICKFILE))
   !      stop
   !   end if
   !   read(LUN,rec=1) surf%mantle
   !   close(LUN)
   !end if

   return
end subroutine io_read_surf
