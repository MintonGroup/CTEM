!**********************************************************************************************************************************
!
!  Unit Name   : regolith_monte_carlo_layer
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Monte carlo method for a strema tube in two layers system  
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_monte_carlo_layer(surfi,erad,deltar,mceb,compmc)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_monte_carlo_layer
   implicit none

   ! Arguments
   type(surftype),intent(in) :: surfi
   real(DP),intent(in)  :: erad,deltar
   real(DP),intent(out) :: mceb,compmc

   ! Traversing a linked list 
   real(DP),parameter :: a = 0.936457 
   real(DP),parameter :: b = 1.12368
   integer(I4B),parameter :: n = 200000000
   real(DP),dimension(3) :: rn,mc_pos
   integer(I4B)       :: i,tot,totlayer
   real(DP)           :: xmax,xmin,ymin,ymax,zmax,zmin,eradi,erado
   real(DP)           :: dmc,costheta,zmc,rmc,vmc,z
   real(DP)           :: totmare,tots

   z = 6.0
   ! Take a stream tube, since it is symmetrical, the streamline running through the body's center is treated as x-axis,
   ! then the width along the center's streamline is varying with a function R(r), "0.5 * rmax * a * tan(b/erad * r)", 
   ! this is y-axis, and the position of a streamtube, depending on the Maxwell Z-model is z-axis. 
   ! Our monte carlo box is simulate the half of the streamtube. So the thickness of a landing ejecta by Monte Carlo method
   ! is twice of the calculated value, but the composition of mare is no change, because the monte carlo mare is as half as the 
   ! the total monte carlo strematube volume.
   ! Determine x range of monte carlo box
   xmax = erad - deltar
   xmin = 0._DP
   ! Determine y range of monte carlo box
   ymax = deltar
   ymin = 0._DP
   ! Determine z range of monte carlo box
   eradi = erad - deltar
   erado = erad + deltar
   zmax = erado / 4.0_DP
   zmin = 0._DP

   tot = 0
   totlayer = 0
   compmc = 0._DP
   mceb = 0._DP

   do i=1,n
      call random_number(rn)
      mc_pos(1) = (xmax - xmin)*rn(1)
      mc_pos(2) = (ymax - ymin)*rn(2)
      mc_pos(3) = (zmax - zmin)*rn(3)

      dmc = 0.5 * deltar * a * tan(b/eradi * mc_pos(1))
      costheta = regolith_quartic_func(mc_pos(1),erad)
      zmc = erad * (1.0 - costheta) * costheta
 
      rmc = sqrt(mc_pos(2)**2 + (mc_pos(3)-zmc)**2)
      
      if (rmc <= dmc) then
           
         tot = tot + 1
                 
         if (mc_pos(3) <= z) then
            totlayer = totlayer + 1
         end if

      end if

   end do

   ! The above calculation does not include the head part of a strematube.
   ! Since the head part of a streamtube intersected with layers can be analytical derived, we can put the volume that is 
   ! calculated from the formula back to the monte carlo result, we should be fine to do get a real ratio of mare that comes 
   ! from a layer in a streamtube. 

   vmc = (xmax-xmin) * (ymax-ymin) * (zmax-zmin) 
   !mceb = 2.0 * real(tot)/real(n) * vmc
   totmare = 0._DP
   tots = 0.0_DP
   call regolith_streamtube_head(surfi,deltar,totmare,tots)
   !compmc = real(totlayer)/real(tot) !without head
   mceb = 2.0 * real(tot)/real(n) * vmc + tots 
   compmc = (2.0*real(totlayer)/real(n) * vmc + totmare)/mceb

   return

end subroutine regolith_monte_carlo_layer
