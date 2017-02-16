!**********************************************************************************************************************************
!
!  Unit Name   : crater_slope_collapse
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Collapses slopes above the angle of repose
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
subroutine crater_slope_collapse(user,surf,crater,domain,deltaMtot)
   use module_globals
   use module_util
   use module_io
   use module_crater, EXCEPT_THIS_ONE => crater_slope_collapse
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(inout) :: deltaMtot
   
   ! Internal variables
   real(DP) :: diffmax,difflim
   real(DP) :: slpmax,slpsq,oldslpmax
   real(DP) :: elchgmax,critical
   real(DP) :: xslp1,xslp2,yslp1,yslp2
   real(DP) :: tslp1sq,tslp2sq,tslp3sq,tslp4sq
   real(DP) :: xgrad,ygrad,tgrad,gradmax
   integer(I4B) :: iradsq
   integer(I4B) :: i,j,inc,incsq,mx1,mx2,mx3,my1,my2,my3,xpi,ypi
   integer(I4B) :: loopcnt,looplim,mult
   real(DP),dimension(:,:),allocatable :: elevarray
   real(DP),dimension(:,:),allocatable :: cumulative_elchange,kappat
   integer(I4B),dimension(:,:,:),allocatable :: indarray
   logical  :: failflag
   character(len=MESSAGESIZE) :: message  ! message for the progress bar
   real(DP) :: rfac

   ! testing
   character(STRMAX) :: filename
   integer(I4B) :: ioerr

   ! Executable Code

   ! Some preliminary setup
   diffmax = 0.25_DP * user%pix**2
   looplim = 10 * crater%fcratpx
   critical = (CRITSLP * user%pix)**2

   !     determine area to effect
   inc = max(min(crater%fcratpx,ceiling(SQRT2*user%gridsize)),1) + 1
   crater%maxinc = max(crater%maxinc,inc)
   incsq = inc**2

   !if (inc >= user%gridsize/2) then
   !   write(*,*) 'D =',crater%fcrat, ' Cl/S =',(crater%fcratpx*user%pix)/domain%side
   !endif
   if (inc >= user%gridsize/2) then
      if (user%testflag) then
          write(*,*) 'Big slope collapse: fcrat =',crater%fcrat, ' Cj/S =',(crater%fcratpx*user%pix)/domain%side
       else
         write(message,'("SlC: Dc=",ES9.2," Cl/S=",F0.3)') crater%fcrat,(crater%fcratpx*user%pix)/domain%side
         call io_updatePbar(message)
       end if
   endif

   allocate(elevarray(-inc:inc,-inc:inc))
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   allocate(indarray(6,-inc:inc,-inc:inc))
   allocate(kappat(-inc:inc,-inc:inc))
   mult = (1 + (inc/(user%gridsize/2)))**2

   cumulative_elchange = 0._DP

   !  begin downslope motion loop
   oldslpmax = huge(oldslpmax)
   kappat = 0.0_DP
   do loopcnt=1,looplim
      failflag = .false.

      ! loop over affected matrix area
      !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
      !$OMP SHARED(inc,loopcnt,incsq,failflag,elevarray,indarray,kappat,diffmax,critical) &
      !$OMP SHARED(crater,user,surf) 
      do j=-inc,inc 
         do i=-inc,inc
            if (loopcnt==1) then
               mx1 = crater%xlpx + i
               my1 = crater%ylpx + j
               mx2 = mx1 + 1
               my2 = my1 + 1
               mx3 = mx1 - 1
               my3 = my1 - 1

               !     periodic boundary conditions 
               call util_periodic(mx1,my1,user%gridsize)
               call util_periodic(mx2,my2,user%gridsize)
               call util_periodic(mx3,my3,user%gridsize)
               indarray(1,i,j) = mx1
               indarray(2,i,j) = my1
               indarray(3,i,j) = mx2
               indarray(4,i,j) = my2
               indarray(5,i,j) = mx3
               indarray(6,i,j) = my3

            else
               mx1 = indarray(1,i,j)
               my1 = indarray(2,i,j)
               mx2 = indarray(3,i,j)
               my2 = indarray(4,i,j)
               mx3 = indarray(5,i,j)
               my3 = indarray(6,i,j)
            end if
            iradsq = i**2 + j**2
            if (iradsq <= incsq) then ! round corners
               ! define location & distance from crater center in pixel space
        
               !     find slopes and gradient
               xslp1 = (surf(mx2,my1)%dem - surf(mx1,my1)%dem) 
               xslp2 = (surf(mx1,my1)%dem - surf(mx3,my1)%dem)
               yslp1 = (surf(mx1,my2)%dem - surf(mx1,my1)%dem)
               yslp2 = (surf(mx1,my1)%dem - surf(mx1,my3)%dem)

               !     find maximum slope
               tslp1sq = (xslp1*xslp1 + yslp1*yslp1)
               tslp2sq = (xslp1*xslp1 + yslp2*yslp2)
               tslp3sq = (xslp2*xslp2 + yslp1*yslp1)
               tslp4sq = (xslp2*xslp2 + yslp2*yslp2)
               slpsq = max(tslp1sq,tslp2sq,tslp3sq,tslp4sq)
               if (slpsq > critical) then
                  kappat(i,j) = diffmax
                  failflag = .true.
               else
                  kappat(i,j) = 0.0_DP
               end if
            end if  
         end do
      end do !     end area loopover
      !$OMP END PARALLEL DO

      if (.not.failflag) exit
     
      elevarray = 0._DP
      
      call util_diffusion_solver(user,surf,2 * inc + 1,indarray(1:2,:,:),kappat,elevarray,mult)

      ! Add the elevation changes back to the dem
      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j)
            ypi = indarray(2,i,j)
            surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + elevarray(i,j)
            surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + elevarray(i,j),0.0_DP)
         end do
      end do

      cumulative_elchange = cumulative_elchange + elevarray

      !***************************************************
      ! TESTING -- makes frames for a slope collapse movie
      !***************************************************
      !if (user%testflag) then
      !   if (mod(loopcnt,1)==0) then
      !      write(filename,'("frames/testprofile_",I4.4,".dat")') loopcnt
      !      open(unit=7,file=filename,status='replace',iostat=ioerr)
      !      if (ioerr==0) then
      !         mx1=crater%xlpx-inc
      !         mx2=crater%xlpx+inc
      !         do i=1,user%gridsize
      !            j=0
      !            if (i>mx1.and.i<mx2) j=1
      !            write(7,*) i*user%pix,surf(i,user%gridsize/2)%dem,j
      !         end do
      !         close(7)
      !      end if
      !   end if
      !end if
      !***************************************************
      

      !     if stuck, get out
      elchgmax = maxval(elevarray)
      if (failflag.and.(abs(elchgmax) < domain%small)) exit

      if (.not.failflag) exit
   end do  ! end downslope motion loops
   
   deltaMtot = deltaMtot + sum(cumulative_elchange)
   deallocate(elevarray,cumulative_elchange,indarray)
   return
   end subroutine crater_slope_collapse
