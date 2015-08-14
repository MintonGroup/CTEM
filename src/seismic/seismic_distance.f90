!
!
!  Unit Name   : seismic_distance
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the maximum distance for seismic shaking
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
subroutine seismic_distance(user,domain,crater,seisdis,maxhits,firstrun)
   use module_globals
   use module_seismic, EXCEPT_THIS_ONE => seismic_distance
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(domaintype),intent(in) :: domain
   type(cratertype),intent(in) :: crater
   real(DP),intent(out) :: seisdis
   integer(I4B),intent(out) :: maxhits
   logical,intent(inout) :: firstrun
   
   ! Internal variables
   real(DP) :: kdiff,startpnt,endpnt,rngdel,gratio,kdiffcalc
   integer(I4B) :: seisdispx,loop
   integer(I4B),parameter :: maxloop = 10000000
   real(DP),parameter :: tol = 1e-12_DP

   ! Executable Code
   kdiff = domain%small
   seisdis = seismic_kdiff_func(user,crater,kdiff,gratio,firstrun,invflag=.true.)
   if (gratio <= 1._DP) then
      startpnt = crater%rad
      endpnt = seisdis
      kdiffcalc = 0._DP
      do loop = 1,maxloop
         rngdel = endpnt - startpnt
         if (abs(rngdel / startpnt) <= tol) exit
         seisdis = startpnt + (rngdel * 0.5_DP)
         kdiffcalc = seismic_kdiff_func(user,crater,seisdis,gratio,firstrun,invflag=.false.)
         if (gratio < 0.1_DP) then
            endpnt = seisdis
         else
            startpnt = seisdis
         endif
      end do
      if (kdiffcalc < kdiff) then
         startpnt = crater%rad
         endpnt = seisdis
         do loop = 1,maxloop
            rngdel = endpnt - startpnt
            if (abs(rngdel / startpnt) <= tol) exit
            seisdis = startpnt + (rngdel * 0.5_DP)
            kdiffcalc = seismic_kdiff_func(user,crater,seisdis,gratio,firstrun,invflag=.false.)
            if (kdiffcalc < kdiff) then
               endpnt = seisdis
            else
               startpnt = seisdis
            endif
         end do
      endif 
   end if
   seisdispx = int(seisdis / user%pix) + 1
   maxhits = (1 + (seisdispx / (user%gridsize / 2)))**2

   return
   end subroutine seismic_distance
