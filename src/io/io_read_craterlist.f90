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
!    Arguments : 
! 
!  Notes       :  only reads needed line [WIP]
!
!**********************************************************************************************************************************
subroutine io_read_craterlist(rclist, user, domain)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_read_craterlist
   implicit none

   ! Arguments
   real(DP),dimension(:,:),intent(out) :: rclist
   type(usertype),intent(inout) :: user
   type(domaintype),intent(in) :: domain

   ! Internals
   integer(I4B)  :: i,ierr
   integer(I4B),parameter :: LUN=7
   real(DP) :: dumm

   ! Executable code

   ! Read the next crater from the craterlist

   open(unit=LUN,file=rcfile,status='old',iostat=ierr)
   if (ierr /= 0) then
      write(*,*) "Unable to open file ",trim(rcfile)
      stop
   end if 

   do i=1,domain%rcnum
      read(LUN,*,iostat=ierr) rclist(1:6,i)
      if (ierr/=0) then
         write(*,*) "Unable to read file ",trim(rcfile)
         stop
      end if
   end do

   close(LUN)
end subroutine io_read_craterlist


