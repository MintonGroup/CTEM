!**********************************************************************************************************************************
!
!  Unit Name   : io_read_const
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Read parameters passed between CTEM and the IDL CTEM driver routine
!
!  Input
!    Arguments : 
!              : 
!
!  Output
!    Arguments : seed  : Current random number generator seeds
!              : ncount : Current step count
!              : totalimpacts  : Total number of impacts in the initial time step
!              : nsteps : Number of steps to execute
!              : restart : F = new run, T = restart old run
! 
!  Notes       : 
!
!**********************************************************************************************************************************

subroutine io_read_const(totalimpacts,ncount,curyear,restart,tallycadence,fracdone,masstot,seedarr)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_read_const
   implicit none
   
   ! Arguments
   integer(I8B),intent(out) :: totalimpacts
   integer(I4B),intent(out) :: ncount,tallycadence
   logical,intent(out) :: restart
   real(DP),intent(out) :: curyear,fracdone,masstot
   integer(I4B),dimension(:),intent(out) :: seedarr
   
   ! Internals
   integer(I4B),parameter :: cfile=21
   integer(I4B) :: ioerr,l

   open(unit=cfile,file=DATFILE,status='old')
   read(cfile,*) totalimpacts,ncount,curyear,restart,tallycadence,fracdone,masstot
   do l=1,size(seedarr)
      read(cfile,*,iostat=ioerr) seedarr(l)
      if (ioerr/=0) seedarr(l)=abs(mod(((l-1)*181)*((seedarr(1)-83)*359),104729))
   end do

   close(cfile)

   
   return

end subroutine io_read_const
