!**********************************************************************************************************************************
!
!  Unit Name   : io_write_surf
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
subroutine io_write_surf(user,surf)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_surf
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf

   ! Internals
   integer(I4B) :: i
   integer(kind=8) :: recsize
   integer(I4B), parameter :: LUN=7
   real(DP) :: dtmp
   real(SP) :: stmp
   integer(I2B) :: itmp


   ! Executable code


   ! Write matrix files
   recsize = sizeof(dtmp) * user%gridsize * user%gridsize
   open(LUN,file=DEMFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) surf%dem
   close(LUN)

   open(LUN,file=EJCOVFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) surf%ejcov
   close(LUN)

   open(LUN,file=DIAMFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   do i=1,user%numlayers 
      write(LUN,rec=i) surf%diam(i)
   end do
   close(LUN)

   open(LUN,file=TIMEFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   do i=1,user%numlayers 
      write(LUN,rec=i) surf%timestamp(i)
   end do
   close(LUN)

   recsize = sizeof(stmp) * user%gridsize * user%gridsize
   open(LUN,file=POSFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   do i=1,user%numlayers 
      write(LUN,rec=2*i-1) surf%xl(i)
      write(LUN,rec=2*i) surf%yl(i)
   end do
   close(LUN)

   if (user%doregotrack) call io_write_regotrack(user,surf) 
   
!   if (user%docrustal_thinning) then
!      recsize = sizeof(itmp) * user%gridsize * user%gridsize
!      open(LUN,file=THICKFILE,status='replace',form='unformatted',recl=recsize,access='direct')
!      write(LUN,rec=1) surf%mantle
!      close(LUN)
!   end if


   return
end subroutine io_write_surf
