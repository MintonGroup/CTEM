!**********************************************************************************************************************************
!
!  Unit Name   : util_chi2
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the reduced chi**2 between two arrays
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine util_filecopy(file1,file2) 
use module_globals
use module_util, EXCEPT_THIS_ONE => util_filecopy
implicit none
character(*),intent(in) :: file1,file2
integer(I4B),parameter :: unit1=20
integer(I4B),parameter :: unit2=22
integer(I4B) :: ioerr
character(len=STRMAX) :: line
character(len=8) :: fmt1,fmt2

open(unit=unit1,file=file1,status='old',iostat=ioerr)
if (ioerr/=0) then
   write(*,*) 'Error opening ',trim(adjustl(file1))
end if
open(unit=unit2,file=file2,status='replace',iostat=ioerr)
if (ioerr/=0) then
   write(*,*) 'Error opening ',trim(adjustl(file2))
end if

write(fmt1,'("(A",I0,")")') STRMAX
do
   read(unit1,fmt=fmt1,iostat=ioerr) line
   if (ioerr/=0) exit
   write(fmt2,'("(A",I0,")")') len_trim(line)
   write(unit2,fmt=fmt2) line
end do


end subroutine util_filecopy
