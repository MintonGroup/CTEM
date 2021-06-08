!**********************************************************************************************************************************
!
!  Unit Name   : io_write_const
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Write parameters passed between CTEM and the IDL CTEM driver routine
!
!  Input 
!    Arguments : 
!              : 
!              : 
!              : 
! 
!  Notes       : 
!
!**********************************************************************************************************************************

subroutine io_write_const(totalimpacts,ncount,curyear,restart,fracdone,masstot,seedarr)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_const
   implicit none
   
   ! Arguments
   integer(I8B),intent(in) :: totalimpacts
   integer(I4B),intent(in) :: ncount
   logical,intent(in) :: restart
   real(DP),intent(in) :: curyear,fracdone,masstot
   integer(I4B),dimension(:),intent(in) :: seedarr
   
   ! Internals
   integer(I4B),parameter :: cfile=21
   integer(I4B) :: l

   open(unit=cfile,file=DATFILE,status='replace')

   1000 format (I17,1X,I12,1X,ES19.12,1X,L1,1X,F9.6,1X,ES19.12)
   write(cfile,1000) totalimpacts,ncount,curyear,restart,fracdone,masstot
   do l=1,size(seedarr)
      write(cfile,'(I12)') seedarr(l)
   end do
   close(cfile)

   return

end subroutine io_write_const
