  !**********************************************************************************************************************************
!
!  Unit Name   : io_write_tally
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes out the tallied crater counts
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


subroutine io_write_tally(tdist,tlist,odist,olist,oposlist,depthdiam,degradation_state)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_tally
   implicit none

   ! Arguments
   real(DP),dimension(:,:),intent(in) :: tdist,tlist,odist
   real(DP),dimension(:),intent(in) :: olist
   real(SP),dimension(:,:),intent(in) :: oposlist
   real(SP),dimension(:),intent(in) :: depthdiam
   real(DP),dimension(:),intent(in) :: degradation_state

   ! Internals
   integer(I4B)  :: i,distl,distc,ioerr,onum
   integer(I4B),parameter :: LUN=7
   logical :: file_exists
   real(DP),dimension(:,:),allocatable :: oldtdist

   ! Executable code

   onum = size(olist)
   distc = size(tdist,1)
   distl = size(tdist,2)
   allocate(oldtdist(distc,distl))
1000 format(3F16.4,2F15.0,F15.6)
2000 format(ES23.15,1X,8(ES14.6,1X,:))
   oldtdist = 0._DP

   inquire(file=tdistfile, exist=file_exists)
   if (file_exists) then
      open(LUN, FILE=TDISTFILE, status='old')
      read(LUN,*,iostat=ioerr) 
      do i=1,distl-1
         read(LUN,*,iostat=ioerr) oldtdist(:,i)
         if (ioerr /= 0) then
            oldtdist(:,i:distl-1) = 0
            exit
         end if
      end do
      close(LUN)
   end if
   open(LUN, file=TDISTFILE, status='REPLACE')
   write(LUN,'("#        Dlo(m)          Dhi(m)        Dmean(m)            dN            N>D               R")')
   do i=1,distl-1
      write(LUN,1000) tdist(1:3,i),oldtdist(4:6,i)+tdist(4:6,i)
   end do
   close(LUN)
   deallocate(oldtdist)

   open(LUN, FILE=ODISTFILE, status='REPLACE')
   write(LUN,'("#        Dlo(m)          Dhi(m)        Dmean(m)            dN            N>D               R")')
   do i=1,distl-1
      write(LUN,1000) odist(:,i)
   enddo
   close(LUN)

   open(LUN, FILE=OLISTFILE, status='REPLACE')
   write(LUN,'("#Dcrat(m)                 xpos(m)        ypos(m)        time(y)        depth/diam        deg_state(m^2)")' )
   do i=1,onum

      write(LUN,2000) olist(i),oposlist(:,i),depthdiam(i),degradation_state(i)
   end do
   close(LUN)


   open(LUN, FILE=TLISTFILE, status='REPLACE')
   write(LUN,'("#Dcrat(m)                 Dimp(m)        xpos(m)        ypos(m)        vimp(m/s)      sinang         time(y)")')
   do i=1,size(tlist,2)
      write(LUN,2000) tlist(:,i)
   end do
   close(LUN)

   return
   end subroutine io_write_tally
