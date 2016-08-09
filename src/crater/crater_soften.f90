!**********************************************************************************************************************************
!
!  Unit Name   : crater_soften
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Softens the terrain of the crater based on empirical studies of equilibrium (see Minton & Fassett 2016)
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
subroutine crater_soften(user,surf,crater,domain)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_soften
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain

   ! Internal variables
   real(DP),dimension(:,:),allocatable :: cumulative_elchange,kappat
   integer(I4B),dimension(:,:,:),allocatable :: indarray
   integer(I4B),parameter :: MAXHITS = 1

   integer(I4B) :: inc,incsq,N,xpi,ypi,iradsq,i,j
   real(DP) :: xp,yp,fradsq,xbar,ybar,areafrac,kappatmax

   kappatmax = SOFTEN_FACTOR * crater%fcrat**SOFTEN_SLOPE
   
   inc = max(min(int(crater%fradpx + 2),PBCLIM * user%gridsize),2) 
   crater%maxinc = max(crater%maxinc,inc)
   incsq = inc**2

   allocate(kappat(-inc:inc,-inc:inc))
   allocate(indarray(2,-inc:inc,-inc:inc))
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   cumulative_elchange = 0._DP
   indarray = inc
   kappat = 0.0_DP

   ! Loop over affected matrix area
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            ! find elevation and grid point
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
            
            xbar = xp - crater%xl 
            ybar = yp - crater%yl

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            
            indarray(1,i,j) = xpi
            indarray(2,i,j) = ypi

            ! We set the diffusion constant to be proportional to the fraction of pixel area covered by the interior of the crater
            areafrac = util_area_intersection(crater%frad,xbar,ybar,user%pix)
             
            kappat(i,j) = kappatmax * areafrac  ! This is the extra per-crater diffusion required to match equilibrium
         end if
      end do
   end do !end area loopover 
   !read(*,*)


   call util_diffusion_solver(user,surf,2 * inc + 1,indarray,kappat,cumulative_elchange,maxhits)
   do j = -inc,inc
      do i = -inc,inc
         xpi = indarray(1,i,j)
         ypi = indarray(2,i,j)
         surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cumulative_elchange(i,j)
         surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cumulative_elchange(i,j),0.0_DP)
      end do
   end do
   deallocate(kappat,indarray,cumulative_elchange)
return
end subroutine crater_soften

