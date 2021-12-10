!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Determines the current degradation state of the crater
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


function crater_get_degradation_state(user,surf,crater,dd) result(Kval)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_get_degradation_state
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   type(cratertype),intent(in) :: crater
   real(SP),intent(out) :: dd
   ! Return value
   real(DP) :: Kval

   ! Internal variables
   integer(I4B) :: inc,xpi,ypi,i,j
   integer(I4B) :: nrim,nbowl,nouter
   real(DP) :: rim,bowl,outer,rad,baseline,ddinit

   ! Definitions of crater dimensions
   real(DP),parameter :: OCUTOFF = 5.5e-2_DP
   real(DP),parameter :: RIMDI = 1.0_DP
   real(DP),parameter :: RIMDO = 1.2_DP
   real(DP),parameter :: BOWLD = 0.2_DP
   real(DP),parameter :: OUTERD = 2.0_DP

   ! Fit for complex crater diffusion model
   real(DP),parameter :: Bconst = 8.0_DP
   real(DP),parameter :: dDpikeC = 1.725_DP ! CTEM version of the Pike (1977) depth-to-diam constant
   real(DP),parameter :: ddsimple = 0.21362307923564000_DP ! CTEM initial depth-to-diameter of simple craterws

   ddinit = min(ddsimple, dDpikeC * (crater%fcrat * 1e-3_DP)**(0.301_DP - 1._DP))

   inc = int(OUTERD * crater%frad / user%pix) + 1 


   nrim = 0
   nbowl = 0
   nouter = 0
   bowl = 0.0_DP
   rim = 0.0_DP
   outer = 0.0_DP

   do j = -inc, inc
      do i = -inc, inc
         rad = sqrt((i**2 + j**2)*1._DP) * user%pix / crater%frad
         baseline = crater%melev + ((i * crater%xslp) + (j * crater%yslp)) * user%pix 
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j
         call util_periodic(xpi,ypi,user%gridsize)
         if (rad <= BOWLD) then
            bowl = bowl + surf(xpi,ypi)%dem - baseline
            nbowl = nbowl + 1 
         else if ((rad >= RIMDI).and.(rad < RIMDO)) then
            rim = rim + surf(xpi,ypi)%dem - baseline
            nrim = nrim + 1 
         else if (rad > RIMDO) then
            outer = outer + surf(xpi,ypi)%dem - baseline
            nouter = nouter + 1 
         end if
      end do
   end do
   rim = rim / nrim 
   bowl = bowl / nbowl 
   outer = outer / nouter
   dd = max((rim - bowl) / crater%fcrat,SMALLFAC)
   Kval = - (1._DP / Bconst * log(real(dd,kind=DP) / ddinit)) * crater%frad**2
   return

end function crater_get_degradation_state

