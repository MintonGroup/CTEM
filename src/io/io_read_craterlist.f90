!**********************************************************************************************************************************
!
!  Unit Name   : io_read_craterlist
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in list of 'real' craters for quasi-MC runs
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : rclist   : List of 'real' craters to read in
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine io_read_craterlist(user, domain)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_read_craterlist
   implicit none

   ! Arguments
   type(usertype),intent(inout) :: user
   type(domaintype),intent(in) :: domain

   ! Internals
   integer(I4B)  :: i,ierr
   integer(I4B),parameter :: LUN=7

   ! Executable code

   ! Read the next crater from the craterlist

   open(unit=LUN,file=rcfile,status='old',iostat=ierr)
   if (ierr /= 0) then
      write(*,*) "Unable to open file ",trim(rcfile)
      stop
   end if 


   do i=1,domain%rcnum
      read(LUN,*,iostat=ierr) user%testimp, user%testvel, user%testang, user%testxoffset, user%testyoffset, user%rctime
      if (ierr/=0) then
         write(*,*) "Unable to read file ",trim(rcfile)
         stop
      end if
   end do
end subroutine io_read_craterlist


