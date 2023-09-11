!**********************************************************************************************************************************
!
!  Unit Name   : io_write_pindex_map
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : write files about pindex map (scoring old melts)
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
subroutine io_write_pindex_map(user,pindex,icrater,ncrat)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_pindex_map
   implicit none

   ! Arguments
   type(usertype),intent(in)                 :: user
   integer(I2B),dimension(:,:),intent(inout) :: pindex
   integer(I8B),intent(in)                   :: icrater
   integer(I8B),intent(in)                   :: ncrat

   ! Regotrack Internals
   integer(I4B), parameter :: LUN=7
   integer(kind=8)         :: recsize
   integer(I2B)            :: itmp
   ! Output multiple "pindex" files
   character(len=255) :: fname
   integer(I2B)       :: n_age
 
   ! Executable code

   n_age   = int(icrater/ncrat)
   recsize = storage_size(itmp) * user%gridsize * user%gridsize / 8
   write(fname,'(a,i6.6)') 'pindex', n_age
   open(LUN,file=fname,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) pindex
   close(LUN)

   return
end subroutine io_write_pindex_map
