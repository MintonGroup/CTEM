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
subroutine crater_slope_collapse(user,surf,crater,domain)
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
   
   ! Internal variables
   real(DP) :: diffmax,difflim
   real(DP) :: slpmax,slpsq,oldslpmax
   real(DP) :: elchange,elchgmax,critical
   real(DP) :: xslp1,xslp2,yslp1,yslp2
   real(DP) :: tslp1sq,tslp2sq,tslp3sq,tslp4sq
   real(DP) :: xgrad,ygrad,tgrad,gradmax
   integer(I4B) :: iradsq
   integer(I4B) :: i,j,inc,incsq,mx1,mx2,mx3,my1,my2,my3
   integer(I4B) :: loopcnt,looplim,mult
   real(DP),dimension(:,:),allocatable :: elevarray
   real(DP),dimension(:,:),allocatable :: gradarray,slparray
   real(DP),dimension(:,:),allocatable :: cumulative_elchange
   integer(I4B),dimension(:,:,:),allocatable :: indarray
   logical  :: failflag,oldfail
   character(len=MESSAGESIZE) :: message  ! message for the progress bar
   real(DP) :: rfac

   ! testing
   !character(STRMAX) :: filename
   !integer(I4B) :: ioerr

   ! Executable Code

   ! Some preliminary setup
   diffmax = 0.25_DP
   looplim = 1000*crater%fcratpx

   !     determine area to effect
   inc = max(min(crater%fcratpx,ceiling(SQRT2*user%gridsize)),1)
   crater%maxinc = max(crater%maxinc,inc)
   incsq = inc*inc

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
   allocate(gradarray(-inc:inc,-inc:inc))
   allocate(slparray(-inc:inc,-inc:inc))
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   allocate(indarray(6,-inc:inc,-inc:inc))

   !mult = inc/(user%gridsize/2) + 1 ! Reduce max diffusion if the loopover area exceeds the periodic boundary conditions
   mult = (1 + (inc/(user%gridsize/2)))**2
   difflim = diffmax/mult

   cumulative_elchange = 0._DP

   !  begin downslope motion loop
   failflag=.true.
   oldslpmax = huge(oldslpmax)
   slparray = 0._DP

   do loopcnt=1,looplim
      !     reset slope failure flag
      !if (mod(loopcnt,100)==0) write(*,*) loopcnt,looplim
      oldfail = failflag
      failflag = .false.

      ! loop over affected matrix area
      !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
      !$OMP SHARED(inc,loopcnt,incsq,difflim,failflag,elevarray,gradarray,slparray,cumulative_elchange,indarray) &
      !$OMP SHARED(crater,user,surf) 
      do j=-inc,inc 
         do i=-inc,inc
            elevarray(i,j) = 0._DP
            gradarray(i,j) = 0._DP
            slparray(i,j) = 0._DP
            if (loopcnt==1) then
               mx1 = crater%xlpx+i
               my1 = crater%ylpx+j
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
            iradsq = i*i+j*j
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
               slparray(i,j) = slpsq
               critical = crater_critical_slope(user,crater,iradsq) * user%pix**2
               if (slpsq > critical) then
                  failflag = .true.
                  ! find elevation change
                  xgrad = (xslp1 - xslp2) 
                  ygrad = (yslp1 - yslp2) 
                  tgrad = xgrad + ygrad
                  elchange = difflim * tgrad

                  ! change digital elevation map
                  elevarray(i,j) = elchange
                  gradarray(i,j) = tgrad
                  cumulative_elchange(i,j) = cumulative_elchange(i,j) + elchange
               end if
            end if  
         end do
      end do !     end area loopover
      !$OMP END PARALLEL DO

     ! See if we have reached maximums
      elchgmax = maxval(elevarray)
      slpmax = sqrt(maxval(slparray))/user%pix
      gradmax = maxval(gradarray)

      ! Add the elevation changes back to the dem
      do j=-inc,inc
         do i=-inc,inc
            mx1=indarray(1,i,j)
            my1=indarray(2,i,j)
            surf(mx1,my1)%dem = surf(mx1,my1)%dem + elevarray(i,j)
            surf(mx1,my1)%ejcov = max(surf(mx1,my1)%ejcov + elevarray(i,j),0.0_DP)
         end do
      end do

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
      if (failflag.and.(abs(elchgmax) < 1.0e-8*domain%small)) exit

      !if (slpmax >= oldslpmax) exit
      !oldslpmax = slpmax

      
      if (.not.failflag) exit
   end do  ! end downslope motion loops
   !j = 0
   !do i = -inc,inc
   !   iradsq = i**2 + j**2
   !   critical = crater_critical_slope(user,crater,iradsq)
   !   write(25,*) i*user%pix/crater%frad,slparray(i,j),critical*user%pix**2
   !end do
   

   deallocate(elevarray,gradarray,slparray,cumulative_elchange,indarray)
   
   return
   end subroutine crater_slope_collapse
