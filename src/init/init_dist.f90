!**********************************************************************************************************************************
!
!  Unit Name   : init_dist
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the distribution arrays
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
subroutine init_dist(user,domain)
   use module_globals
   use module_init, EXCEPT_THIS_ONE => init_dist
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(domaintype),intent(inout) :: domain

   ! Internals
   integer(I4B), parameter :: LUN = 7
   integer(I4B)           :: ierr
   real(DP)               :: testreal

   ! Executable code

   ! Get size of crater production function array
   open(unit=LUN,file=trim(adjustl(user%sfdfile)),status="old",iostat=ierr)
   if (ierr /= 0) then
      write(*,*) "Unable to open file ", trim(user%sfdfile)
      stop
   end if

   domain%pnum = 0
   do
      read(LUN,*,iostat=ierr) testreal,testreal
      if (ierr/=0) exit
      domain%pnum = domain%pnum + 1
   end do

   if (domain%pnum == 0) then
      write(*,*) "No valid entries in ",trim(adjustl(user%sfdfile))
   end if
   close(LUN)

   ! Get size of velocity distribution array
   open(unit=LUN,file=trim(adjustl(user%velfile)),status="old",iostat=ierr)
   if (ierr /= 0) then
   write(*,*) "Unable to open file ",trim(user%velfile)
      stop
   end if

   domain%vnum = 0
   do
      read(LUN,*,iostat=ierr) testreal,testreal,testreal
      if (ierr /= 0) exit
      domain%vnum = domain%vnum + 1
   end do
   if (domain%vnum == 0) then
      write(*,*) "No valid entries in ",trim(adjustl(user%velfile))
   end if
   close(LUN)

   ! Get size of real crater list array
   if (user%doquasimc) then
      open(unit=LUN,file=trim(adjustl(rcfile)),status="old",iostat=ierr)
      if (ierr /= 0) then
         write(*,*) "Unable to open file ", trim(adjustl(rcfile))
         stop
      end if

      domain%rcnum = 0
      do
         read(LUN,*,iostat=ierr) testreal,testreal,testreal,testreal,testreal,testreal
         if (ierr/=0) exit
         domain%rcnum = domain%rcnum + 1
      end do

      if (domain%rcnum == 0) then
         write(*,*) "No valid entries in ",trim(rcfile)
      end if
      close(LUN)
   else
      domain%rcnum = 1
   end if


   return
end subroutine init_dist
