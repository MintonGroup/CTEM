  !**********************************************************************************************************************************
!
!  Unit Name   : io_write_dist
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes the impactor production and scaled crater arrays
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       :  
!
!**********************************************************************************************************************************


subroutine io_write_dist(pdist,crtscl,domain,mass)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_dist
   implicit none

   ! Arguments
   real(DP),dimension(:,:),intent(in) :: pdist,crtscl
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: mass

   ! Internals
   integer(I4B)  :: i
   integer(I4B),parameter :: LUN=7

   ! Executable code


   open(LUN, FILE=PDISTFILE, status='UNKNOWN')
   write(LUN,'("#        Dlo(m)          Dhi(m)        Dmean(m)            dN            N>D               R")')
   do i=1,domain%pdistl
      write(LUN,'(6(ES17.10,1X))') pdist(:,i)
   enddo
   close(LUN)

   open(LUN, FILE=CRTSCLFILE, status='UNKNOWN')
   write(LUN,'("#Dimp(m)           Dcrat(m)")')
   do i=1,domain%pnum
      write(LUN,'(2(ES17.10))') crtscl(:,i)
   enddo
   close(LUN)

   open(LUN, FILE=MASSFILE, status='UNKNOWN')
   write(LUN,*) mass
   close(LUN)

   return

   end subroutine io_write_dist
