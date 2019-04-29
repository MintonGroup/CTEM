!**********************************************************************************************************************************
!
!  Unit Name   : io_updatePbar
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : updates the progress bar
!
!  Input
!    Arguments : infile : input filename
!
!  Output
!    Arguments : 
! 
!  Notes       : 
!
!**********************************************************************************************************************************
subroutine io_updatePbar(message)
use module_globals
use module_io, EXCEPT_THIS_ONE => io_updatePbar
implicit none
character(len=*),intent(in) :: message
integer(I4B) ::k,endpoint
integer(I4B) :: perc
character(len=PBARSIZE+7) ::bar
integer(I4B),save :: flip
character(len=STRMAX) :: fmtlabel 
character(len=MESSAGESIZE),save :: persistent_message


bar="???% ("
bar(7:6+PBARSIZE)=" "
bar(7+PBARSIZE:7+PBARSIZE)= ")"
if (message .ne. "") then
   persistent_message = message
end if

perc=floor(100*real(pbarpos)/real(PBARRES))
write(unit=bar(1:3),fmt="(i3)") perc
endpoint = min(ceiling(real(pbarpos*PBARSIZE)/real(PBARRES)),PBARSIZE)
do k = 1,endpoint !- 1 
   bar(6+k:6+k)=pbarchar(k:k)
end do
!select case(flip)
!case(1)
!   bar(6+endpoint:6+endpoint)="/"
!case(2)
!   bar(6+endpoint:6+endpoint)="-"
!case(3)
!   bar(6+endpoint:6+endpoint)="\"
!case(4)
!   bar(6+endpoint:6+endpoint)="|"
!end select
!flip = flip + 1
if (flip > 4) flip = 1
! print the progress bar.  
write(fmtlabel,'("(A1,A",I2.2,",1X,A",I2.2,",$)")') PBARSIZE+7,MESSAGESIZE
write(*,fmt=fmtlabel) char(13), bar,persistent_message !trim(adjustl(persistent_message))

return
end subroutine io_updatePbar

