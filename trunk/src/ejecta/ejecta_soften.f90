!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_soften
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Softens the terrain under the ejecta  using a box filter model where the size of the box is proportional to the 
!                thickness of the ejecta  
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
subroutine ejecta_soften(user,surf,N,indarray,cumulative_elchange)
   use module_globals
   use module_util
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_soften
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   integer(I4B),intent(in) :: N
   integer(I4B),dimension(2,N,N),intent(in) :: indarray
   real(DP),dimension(N,N),intent(inout) :: cumulative_elchange 
   !real(DP),parameter :: SOFTEN_FACTOR = 2.50_DP ! Constant in topographic diffusion term for ejecta blanket softening
   real(DP),parameter :: SOFTEN_FACTOR = 2.00_DP ! Constant in topographic diffusion term for ejecta blanket softening

   ! Internal variables
   integer(I4B) :: maxhits
   real(DP),dimension(N,N) :: ebharr
   real(DP),dimension(N,N) :: kdiff


   ! Save the original ejecta blanket thickness
   ebharr = cumulative_elchange 

   maxhits = 1

   ! Diffusion constant for 1 time unit was found to be proportional to ejecta thickness times the pixel size
   kdiff = SOFTEN_FACTOR * user%pix * ebharr 
   !TESTING
   !kdiff = user%sf * user%pix**user%k * ebharr**user%p !**0.900

   call util_diffusion_solver(user,surf,N,indarray,kdiff,cumulative_elchange,maxhits)

   ! Add back the overall shape of the ejecta
   cumulative_elchange = cumulative_elchange + ebharr

return
end subroutine ejecta_soften

