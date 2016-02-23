!**********************************************************************************************************************************
!
!  Unit Name   : init_regolith_parab
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the shape of parabola in the regolith stack
!  
!
!  Input
!    Arguments : regolayer
!
!  Output
!    Arguments : regolayer
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine init_regolith_parab(user,surf)
   use module_globals
   use module_regolith
   use module_util
   use module_init, EXCEPT_THIS_ONE => init_regolith_parab
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(regolayertype) :: mare,highland,porous,finalc
   integer(I4B) :: i,j,xpi,ypi,maresize

   ! Initial parabola shape parameters
   real(DP), parameter :: TRDDRATIO = (1.0_DP/3.0_DP + 1.0_DP/4.0_DP) / 2.0_DP
   integer(I4B), parameter :: inc = 534
   integer(I4B) :: incsq, xlpx, ylpx,iradsq
   real(DP) :: fradsq, xl, yl, xp, yp, lradsq
   real(DP) :: parab, vcorr, cform, frad !final crater parabola profile
   real(DP) :: trparab, trvcorr, trdepth, trform, rad !transient crater parabola profile
   real(DP) :: porous_thick, finalc_depth, porous_depth

   ! One layer first   
   do j = 1,user%gridsize
      do i = 1,user%gridsize
            highland%thickness = 1000.00_DP
            highland%meltfrac  = 0._DP 
            highland%comp      = 0._DP
            call regolith_push(surf(i,j),highland)
      end do
   end do

   ! Start to build parabola shaped profile on push-pop system
   ! Calculate the thickness of a layer for each pixel between the transient crater
   ! and the final crater. 
   ! Final crater's parabola shape parameters
   vcorr = 172.07029707833178
   parab = 6.6693672265677364E-004
   ! Transient crater's parabola shape parameters
   rad   = 455.81535649889213
   trdepth = TRDDRATIO * rad * 2.0_DP
   trparab = trdepth / (rad**2)
   trvcorr = trdepth
   
   incsq = inc**2
   xlpx  = 250
   ylpx  = 250
   xl    = 5000.00000
   yl    = 5000.00000
   fradsq = 324636.93628158141

   do j=-inc,inc
      do i=-inc,inc
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            xpi = xlpx + i
            ypi = ylpx + j
            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
            lradsq = (xl - xp)**2 + (yl - yp)**2
            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            if (lradsq < fradsq) then
               cform = vcorr - (parab * lradsq)
               trform  = trparab * lradsq - trvcorr
               if (trform<(-cform) .and. trform < 0._DP) then
                  porous_thick = abs(abs(trform) - abs(cform))
                  finalc_depth = abs(cform)
                  porous_depth = finalc_depth + porous_thick
                  call regolith_traverse_pop(-1.0_DP * porous_depth,surf(xpi,ypi))
                  ! Push the porous layer first
                  porous%thickness = porous_thick * 0.8
                  porous%comp      = 1.0_DP
                  porous%meltfrac  = 0.0_DP
                  call regolith_push(surf(xpi,ypi),porous)
                  porous%thickness = porous_thick * 0.2
                  porous%comp      = 0.0_DP
                  porous%meltfrac  = 0.0_DP
                  call regolith_push(surf(xpi,ypi),porous)  
                  finalc%thickness = finalc_depth
                  finalc%comp      = 0.0_DP
                  finalc%meltfrac  = 0.0_DP
                  call regolith_push(surf(xpi,ypi),finalc)
               end if
            end if
          end if
      end do
   end do

   return
end subroutine init_regolith_parab

