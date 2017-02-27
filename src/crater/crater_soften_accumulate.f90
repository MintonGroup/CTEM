!**********************************************************************************************************************************
!
!  Unit Name   : crater_soften_accumulate
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Softens the terrain of the crater based on empirical studies of equilibrium (see Minton & Fassett 2016).
!                Accumulates the diffusivity of the grid for later application in order to improve performance
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
subroutine crater_soften_accumulate(user,surf,crater,domain,kdiff)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_soften_accumulate
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(inout) :: kdiff

   ! Internal variables
   integer(I4B) :: inc,cinc,incsq,N,xpi,ypi,iradsq,i,j
   real(DP) :: kappatmax,lrad,lradsq,xp,yp,fradsq,areafrac,xbar,ybar
   real(DP),dimension(user%gridsize,user%gridsize) :: craterhole,kdifftmp
   logical,dimension(user%gridsize,user%gridsize) :: hit 

   hit = .false.

   kappatmax = user%soften_factor / (PI * user%soften_size**2) * crater%frad**(user%soften_slope - 2.0_DP)
   !Take out the intrinsic crater erosion contribution 
   kappatmax = max(kappatmax - 0.84_DP / (PI * user%soften_size**2) * crater%frad**2,0.0_DP) 

   inc = int(min(user%soften_size * crater%frad / user%pix, user%gridsize / SQRT2)) + 2 
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   incsq = inc**2
   cinc = crater%fradpx + 1
   craterhole = 1.0_DP
   ! Cut out a hole for the new crater so we don't soften it with itself
   do j = -cinc,cinc
      do i = -cinc,cinc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! Find distance from crater center to current pixel center in real space
         xp = xpi * user%pix
         yp = ypi * user%pix

         xbar = xp - crater%xl 
         ybar = yp - crater%yl
   
         areafrac = util_area_intersection(crater%frad,xbar,ybar,user%pix)

         call util_periodic(xpi,ypi,user%gridsize)
         craterhole(xpi,ypi) = 1.0_DP - areafrac
      end do
   end do

   ! Loop over affected matrix area
   !!$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !!$OMP SHARED(user,crater,fradsq,inc,incsq,kappatmax,craterhole,nothit) &
   !!$OMP REDUCTION(+:kdiff)
   do j = -inc,inc  ! Do the loop in pixel space
      do i = -inc,inc
         ! find distance from crater center
         iradsq = i**2 + j**2
         if (iradsq <= incsq) then
            ! find elevation and grid point
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix

            xbar = xp - crater%xl 
            ybar = yp - crater%yl
            
            lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            areafrac = util_area_intersection(user%soften_size * crater%frad,xbar,ybar,user%pix)
            areafrac = areafrac * craterhole(xpi,ypi)
            
            if (.not.hit(xpi,ypi)) then 
               kdiff(xpi,ypi) = kdiff(xpi,ypi) + kappatmax * areafrac 
               hit = .true.
            end if

         end if

      end do
   end do !end area loopover 
   !!$OMP END PARALLEL DO
return
end subroutine crater_soften_accumulate

