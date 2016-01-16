!**********************************************************************************************************************************
!
!  Unit Name   : seismic_table_define
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
subroutine seismic_table_define(user,domain,crater,srim,totdiff,diffmax,srng,kdifflim,diffcnt,firstrun)
   use module_globals
   use module_seismic, EXCEPT_THIS_ONE => seismic_table_define
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(domaintype),intent(in) :: domain
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: srim,diffmax,totdiff
   integer(I4B),intent(out) :: diffcnt
   integer(I4B),dimension(:),allocatable,intent(out) :: srng
   real(DP),dimension(:),allocatable,intent(out) :: kdifflim
   logical,intent(inout) :: firstrun
   
   ! Internal variables
   real(DP) :: kdiff,lrad,startpnt,endpnt,rngdel,gratio,kdiffcalc,difflim
   integer(I4B) :: srngmul,i,loop,seisdis,maxcnt
   integer(I4B),parameter :: maxloop = 10000000
   real(DP),parameter :: tol = 1e-12_DP

   ! Executable Code
   kdiff = domain%small
   do i=1,maxloop
      lrad = seismic_kdiff_func(user,crater,kdiff,gratio,firstrun,invflag=.true.)
      if (gratio <= 1._DP) then
         startpnt = srim
         endpnt = lrad
         kdiffcalc = 0._DP
         do loop=1,maxloop
            rngdel = endpnt - startpnt
            if (abs(rngdel/startpnt) <= tol) exit
            lrad = startpnt + (rngdel * 0.5_DP)
            kdiffcalc = seismic_kdiff_func(user,crater,lrad,gratio,firstrun,invflag=.false.)
            if (gratio < 0.1_DP) then
               endpnt = lrad
            else
               startpnt = lrad
            endif
         end do
         if (kdiffcalc < kdiff) then
            startpnt = srim
            endpnt = lrad
            do loop=1,maxloop
               rngdel = endpnt - startpnt
               if (abs(rngdel / startpnt) <= tol) exit
               lrad = startpnt + (rngdel*0.5_DP)
               kdiffcalc = seismic_kdiff_func(user,crater,lrad,gratio,firstrun,invflag=.false.)
               if (kdiffcalc < kdiff) then
                  endpnt = lrad
               else
                  startpnt = lrad
               endif
            end do
         endif 
         !kdiff = kdiffcalc
      end if
      seisdis = int(lrad/user%pix) + 1
      srngmul = (1 + (seisdis/(user%gridsize/2)))**2
      difflim = diffmax/srngmul
      if (i==1) then
         maxcnt = int(min(totdiff/difflim,real(huge(maxcnt)-1,kind=DP)) )+ 1
         allocate(srng(maxcnt))
         allocate(kdifflim(maxcnt))
      end if
      kdiff = kdiff + difflim
      srng(i) = seisdis
      kdifflim(i) = difflim
      diffcnt = i
      if (kdiff > totdiff) exit
   end do 
   return
   end subroutine seismic_table_define
