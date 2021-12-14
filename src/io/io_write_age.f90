!**********************************************************************************************************************************
!
!  Unit Name   : io_write_age
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for terrain grids for age information
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
subroutine io_write_age(user,surf,n_size,icrater,ncrat)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_age
   implicit none

   ! Arguments
   type(usertype),intent(in)                :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   integer(I2B),intent(in)                  :: n_size
   integer(I8B),intent(in)                  :: icrater
   integer(I8B),intent(in)                  :: ncrat

   ! Regotrack Internals
   integer(I4B) :: i,j,k
   integer(I4B), parameter :: LUN=7
   real(DP),parameter                                  :: sdepth = 0.1
   real(SP),dimension(user%gridsize,user%gridsize)     :: agetop
   real(SP),dimension(:),allocatable                   :: age_prev
   real(SP),dimension(:),allocatable                   :: agedepthtot
   integer(kind=8)                                     :: recsize
   real(SP)                                            :: stmp
   real(SP)                                            :: agetot, age_weighted
   type(regolisttype),pointer                          :: current => null()
   real(DP)                                            :: depth, depth_prev
   real(SP)                                            :: recyclratio
   ! Output multiple "comphisto" files
   character(len=255) :: fname
   integer(I4B)       :: n_age
 
   ! Executable code
   allocate(age_prev(n_size))
   allocate(agedepthtot(n_size))
   agetop = 0.0_SP
   
   do j=1,user%gridsize
      do i=1,user%gridsize
         
         depth       = 0.0 
         depth_prev  = 0.0 
         current => surf(i,j)%regolayer
         age_prev(:) = current%regodata%age(:)
         agedepthtot = 0.0_SP
         do 
          if (depth > sdepth) then
             recyclratio    = real( (depth - sdepth) / (depth - depth_prev) )
             agedepthtot(:) = agedepthtot(:) - recyclratio * age_prev(:)
             depth          = sdepth
             exit
           end if
           depth_prev     = depth
           depth          = depth + current%regodata%thickness
           agedepthtot(:) = agedepthtot(:) + current%regodata%age(:)
           age_prev(:)    = current%regodata%age(:)
           current => current%next
         end do

         age_weighted = 0.0_SP
         do k=1, n_size
            age_weighted = age_weighted + agedepthtot(k) * (float(k) - 0.5) * 0.5 
         end do

         agetot       = sum(agedepthtot(:))
         age_weighted = age_weighted / agetot
         if (agetot == 0.0_SP .and. age_weighted == 0.0_SP) then            
            agetop(i,j) = -0.5
         else 
            agetop(i,j) = age_weighted
         end if 

      end do 
   end do

   deallocate(age_prev,agedepthtot)

   n_age   = max(int(icrater/ncrat), 1)
   recsize = sizeof(stmp) * user%gridsize * user%gridsize
   write(fname,'(a,i6.6)') 'agetop', n_age
   open(LUN,file=fname,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) agetop
   close(LUN)
  
   return
end subroutine io_write_age
