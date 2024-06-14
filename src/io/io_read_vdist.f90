  !**********************************************************************************************************************************
!
!  Unit Name   : io_read_vdist
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in the impactor vdistuction function array
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : vdist   : Impactor vdistuction function
! 
!  Notes       :  
!
!**********************************************************************************************************************************


subroutine io_read_vdist(vdist,user,domain)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_read_vdist
   implicit none

   ! Arguments
   real(DP),dimension(:,:),intent(out) :: vdist
   type(usertype),intent(in) :: user
   type(domaintype),intent(in) :: domain

   ! Internals
   integer(I4B)  :: i,ierr
   integer(I4B),parameter :: LUN=7

   ! Executable code

   ! Read in velocity distribution file
   open(unit=LUN,file=trim(adjustl(user%velfile)),status='old',iostat=ierr)
   if (ierr /= 0) then
      write(*,*) "Unable to open file ",trim(adjustl(user%velfile))
      stop
   end if   

   do i=1,domain%vnum
      read(LUN,*,iostat=ierr) vdist(1:3,i)
      if (ierr/=0) then
         write(*,*) "Unable to read file ",trim(adjustl(user%velfile))
         stop
      end if
   end do

   close(LUN)

   return

   end subroutine io_read_vdist
