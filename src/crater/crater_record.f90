!**********************************************************************************************************************************
!
!  Unit Name   : crater_record
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : records the new crater in an available layer
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine crater_record(user,surf,crater,melev,xslp,yslp)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_record
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   real(DP),intent(in) :: melev,xslp,yslp

   ! Internal variables
   real(DP) :: rimdis,rimdissq,avgrim,baseline
   real(SP) :: depth
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq,nrim,nbowl
   integer(I2B) :: isrim
   real(DP),dimension(:),allocatable :: rimelevation,bowlelevation
   integer(I4B),dimension(:),allocatable :: elevind
   

   ! Executable code

   ! determine area to effect
   rimdis = (crater%frad / user%pix)
   inc  = max(int(rimdis) + 1, 1)
   incsq = inc**2
   allocate(rimelevation(4 * incsq))
   allocate(bowlelevation(4 * incsq))
   allocate(elevind(4 * incsq))

   rimdissq = rimdis**2

   ! Loop over affected matrix area
   nrim = 0
   nbowl = 0
   do j = -inc,inc
      do i = -inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
            baseline = melev + ((i * xslp) + (j * yslp)) * user%pix

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! find out if we're part of the counting rim or not
            if ((iradsq * 1._DP) >= rimdissq) then
               nrim = nrim + 1
               rimelevation(nrim) = surf(xpi,ypi)%dem - baseline
            else
               nbowl = nbowl + 1
               bowlelevation(nbowl) = surf(xpi,ypi)%dem - baseline
            end if
         end if

      end do
   end do !end area loopover 

   if (nrim < 2) then
      write(*,*) 'Crater rim too small to count! Increase the value of SMALLESTCOUNTABLE or rimfrac'
      avgrim = -huge(depth)
   else
      avgrim = sum(rimelevation(1:nrim)) / real(nrim, kind = DP)
   end if

   do j=-inc,inc
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            baseline = melev + ((i * xslp) + (j * yslp)) * user%pix

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! find out if we're part of the counting rim or not
            if ((iradsq*1._DP) >= rimdissq) then
               isrim = 1
            else
               isrim = 0
            end if
            ! Calculate the current depth below the average rim height
            depth = real(avgrim - (surf(xpi,ypi)%dem - baseline), kind = SP)

            ! record crater in available layer
            call util_add_to_layer(user,surf(xpi,ypi),isrim,crater%fcrat,crater%xl,crater%yl,depth,real(baseline, kind = SP))
         end if

      end do
   end do 

   deallocate(rimelevation,bowlelevation,elevind)
   return
end subroutine crater_record

