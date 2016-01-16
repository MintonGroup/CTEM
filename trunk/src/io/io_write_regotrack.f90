!**********************************************************************************************************************************
!
!  Unit Name   : io_write_regotrack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for terrain grids for regolith tracking        
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
subroutine io_write_regotrack(user,surf)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_regotrack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf

   ! Regotrack Internals
   integer(I4B) :: i,j,k
   integer(I4B), parameter :: LUN=7
   integer(I4B), parameter :: LUM=8
   integer(I4B), parameter :: LUC=9
   type(regolayertype),pointer :: current
   real(DP),dimension(user%gridsize,user%gridsize) :: regotop,comp,melt
   integer(I4B),dimension(user%gridsize,user%gridsize) :: stacks_num
   integer(kind=8) :: recsize
   real(DP) :: dtmp
   integer(I4B) :: itmp
   real(DP),dimension(user%gridsize,user%gridsize) :: comptop!,surface
   real(DP),dimension(:),allocatable :: marehisto
   real(DP) :: mare, z

   ! Mixing 
   !real(DP),parameter :: zmix = 0.0_DP 
   !real(DP) :: z, zmare

   ! Executable code
   open(LUN,file=MELTFILE,status='replace',form='unformatted')
   open(LUM,file=REGOFILE,status='replace',form='unformatted')
   open(LUC,file=COMPFILE,status='replace',form='unformatted')

   do j=1,user%gridsize
      do i=1,user%gridsize
         stacks_num(i,j) = 0
         current => surf(i,j)%regolayer
         comptop(i,j) = current%comp
         do 
          if (.not. associated(current%next)) exit
          stacks_num(i,j) = stacks_num(i,j) + 1
          regotop(i,j) = current%thickness
          comp(i,j) = current%comp
          melt(i,j) = current%meltfrac
          write(LUM) regotop(i,j)
          write(LUC) comp(i,j)
          write(LUN) melt(i,j) 
          !if (regotop(i,j) == 0._DP) then
          !write(*,*) i,j,stacks_num(i,j),regotop(i,j)!,comp(i,j)
          !end if
          current => current%next
         end do
      end do 
   end do
   close(LUM)
   close(LUC)

   recsize = sizeof(dtmp) * user%gridsize * user%gridsize
   open(LUN,file='comptop.dat',status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) comptop
   close(LUN)

   allocate(marehisto(user%gridsize))
   open(LUN,file='comphisto',status='replace')
   do i=1,user%gridsize
      marehisto(i) = 0.0_DP
      mare = 0.0_DP
      do j=1,user%gridsize
         mare = mare + comptop(i,j)
      end do
      marehisto(i) = mare/real(user%gridsize)
      write(LUN,*) real(i-user%gridsize/2)*user%pix/1000.0,marehisto(i)*100.0
   end do
   close(LUN)
   deallocate(marehisto)

   recsize = sizeof(itmp) * user%gridsize * user%gridsize
   open(LUN,file=STACKNUMFILE,status='replace',form='unformatted',recl=recsize,access='direct')
   write(LUN,rec=1) stacks_num
   close(LUN)

   return
end subroutine io_write_regotrack
