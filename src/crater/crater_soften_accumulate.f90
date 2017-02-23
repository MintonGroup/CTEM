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
   integer(I4B) :: SOFTEN_SIZE = 100 ! Constant in topographic diffusion term for crater softening

   integer(I4B) :: inc,incsq,N,xpi,ypi,iradsq,i,j
   real(DP) :: kappatmax,lrad,lradsq,xp,yp,fradsq,areafrac,xbar,ybar

   ! TESTING
   !open(unit=55,file="SOFTEN_FACTOR.test",status="old")
   !read(55,*) SOFTEN_FACTOR
   !close(55)
   !********


   kappatmax = user%soften_factor / (PI * (SOFTEN_SIZE * crater%frad)**(user%soften_slope - 2.0_DP))
   kappatmax = kappatmax - 0.84_DP / (PI * (SOFTEN_SIZE * crater%frad)**2.0)
   inc = SOFTEN_SIZE * crater%fradpx + 2 
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   incsq = inc**2

   ! Loop over affected matrix area
   !!$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !!$OMP SHARED(user,crater,fradsq,inc,incsq,kappatmax,SOFTEN_SIZE) &
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
            areafrac = util_area_intersection(SOFTEN_SIZE * crater%frad,xbar,ybar,user%pix)
            
            if (lradsq > fradsq) then 
               kdiff(xpi,ypi) = kdiff(xpi,ypi) + kappatmax * areafrac 
            end if

         end if

      end do
   end do !end area loopover 
   !!$OMP END PARALLEL DO

return
end subroutine crater_soften_accumulate

