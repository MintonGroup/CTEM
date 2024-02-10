  !**********************************************************************************************************************************
!
!  Unit Name   : io_read_prod
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in the impactor production function array
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : prod   : Impactor production function
! 
!  Notes       :  
!
!**********************************************************************************************************************************


subroutine io_read_prod(prod,user,domain)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_read_prod
   implicit none

   ! Arguments
   real(DP),dimension(:,:),intent(out) :: prod
   type(usertype),intent(in) :: user
   type(domaintype),intent(in) :: domain

   ! Internals
   integer(I4B)  :: i,ierr
   integer(I4B),parameter :: LUN=7

   ! Executable code

   ! Read in the size-frequency distribution file
   open(unit=LUN,file=trim(adjustl(user%sfdfile)),status='old',iostat=ierr)
   if (ierr /= 0) then
      write(*,*) "Unable to open file ",trim(user%sfdfile)
      stop
   end if   

   do i=1,domain%pnum
      read(LUN,*,iostat=ierr) prod(1:2,i)
      if (ierr/=0) then
         write(*,*) "Unable to read file ",trim(user%sfdfile)
         stop
      end if
   end do

   close(LUN)

   return

   end subroutine io_read_prod
