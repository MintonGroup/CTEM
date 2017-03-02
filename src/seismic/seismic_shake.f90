!**********************************************************************************************************************************
!
!  Unit Name   : seismic_shake
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Implements seismic shaking for a given crater
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
subroutine seismic_shake(user,surf,crater,domain)
   use module_globals
   use module_crater
   use module_util
   use module_io
   use module_seismic, EXCEPT_THIS_ONE => seismic_shake
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   
   ! Internal variables
   real(DP) :: totdiff,lradsq,radsq,gratio,xp,yp,seisdis
   integer(I4B) :: iradsq
   integer(I4B) :: i,j,inc,incsq,xpi,ypi
   integer(I4B) :: maxhits
   real(DP),dimension(:,:),allocatable :: cumulative_elchange
   real(DP),dimension(:,:),allocatable :: kdiff
   integer(I4B),dimension(:,:,:),allocatable :: indarray
   logical :: firstrun
   character(len=MESSAGESIZE) :: message  ! message for the progress bar

   ! Executable Code

   ! Some preliminary setup

   ! find seismic constant terms for this crater
      crater%kdiffterm = SHEFF * (user%seisq**QFAC) * (user%neff**NFAC) * (crater%imp**PFAC)
      crater%kdiffterm = crater%kdiffterm * (crater%impvel**VFAC) * (user%gaccel**GFAC)
      crater%saccelterm = user%neff * user%prho * (crater%impvel**2) * (SEISFREQ**2) * (crater%imp**3)

   !  find acceleration ratio at the crater rim
   totdiff = seismic_kdiff_func(user,crater,crater%frad,gratio,invflag=.false.)

   ! apply seismic diffusion if magnitude is great enough
   if (totdiff < domain%small) return

   call seismic_distance(user,domain,crater,seisdis,maxhits)

   !     report to screen
   if (seisdis / user%pix > (user%gridsize/2)) then
      if (user%testflag) then
          write(*,*) 'Big Seismic: fcrat =',crater%fcrat, ' Se/S =',seisdis / domain%side, ' Serim =',totdiff
       else
         write(message,'("Seis: Dc=",ES9.2," Se/S=",F0.2)') crater%fcrat,seisdis / domain%side
         call io_updatePbar(message)
       end if
   endif
 
   !     determine area to effect
   inc = max(min(int(seisdis / user%pix) + 1,PBCLIM*user%gridsize),1)

   crater%maxinc = max(crater%maxinc,inc)
   incsq = inc**2
   radsq = crater%rad**2

   allocate(indarray(2,-inc:inc,-inc:inc))
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   allocate(kdiff(-inc:inc,-inc:inc))

   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !$OMP SHARED(user,crater,inc,incsq,indarray,radsq,kdiff) 
   do j = -inc,inc 
      do i = -inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! Find distance from crater center to current pixel center in real space
         xp = xpi * user%pix
         yp = ypi * user%pix

         lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)

         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi

         if ((iradsq <= incsq) .and. (lradsq >= radsq)) then
            kdiff(i,j) = seismic_kdiff_func(user,crater,sqrt(lradsq),gratio,invflag=.false.)
         else if (lradsq < radsq) then
            kdiff(i,j) = seismic_kdiff_func(user,crater,sqrt(radsq),gratio,invflag=.false.)
         else
            kdiff(i,j) = 0.0_DP
         end if
      end do
   end do
   !$OMP END PARALLEL DO

   call util_diffusion_solver(user,surf,2 * inc + 1,indarray,kdiff,cumulative_elchange,maxhits)
    
   ! Add the total diffusion to the layers
   do j=-inc,inc
      do i=-inc,inc
         xpi = indarray(1,i,j)
         ypi = indarray(2,i,j)
         surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cumulative_elchange(i,j)
         surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cumulative_elchange(i,j),0.0_DP)
      end do
   end do

   deallocate(cumulative_elchange,indarray,kdiff)
   return

   end subroutine seismic_shake
