!**********************************************************************************************************************************
!
!  Unit Name   : io_crater_profile
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for terrain grids           
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : surf : Surface expression matrix
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine io_crater_profile(user,surf)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_crater_profile
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf

   ! Internals
   integer(I4B) :: i,j
   integer(I4B), parameter :: LUN=7

   ! Executable code
   write(*,*) 'Writing crater profile to testprofile.dat'
   j=size(surf,2)/2
   open(unit=LUN,file='testprofile.dat',status='replace')
   do i=1,size(surf,1)
      write(LUN,*) user%pix*i,surf(i,j)%dem
   end do

   close(LUN)

   return
end subroutine io_crater_profile
