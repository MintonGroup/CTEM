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
subroutine crater_soften(user,surf,crater,domain,kdiff)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_soften
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(in) :: kdiff

   ! Internal variables
   real(DP),dimension(0:user%gridsize + 1,0:user%gridsize + 1) :: cumulative_elchange
   integer(I4B),dimension(2,0:user%gridsize + 1,0:user%gridsize + 1) :: indarray
   integer(I4B) :: i,j,xpi,ypi
   integer(I4B),parameter :: MAXHITS = 1

   ! Create box for soften calculation (will be no bigger than the grid itself)
   do j = 0,user%gridsize + 1
      do i = 0,user%gridsize + 1
         xpi = i
         ypi = j
         call util_periodic(xpi,ypi,user%gridsize)
         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi
      end do
   end do
   call util_diffusion_solver(user,surf,user%gridsize + 2,indarray,kdiff,cumulative_elchange,maxhits)
   do j = 1,user%gridsize
      do i = 1,user%gridsize
         surf(i,j)%dem = surf(i,j)%dem + cumulative_elchange(i,j)
         surf(i,j)%ejcov = max(surf(i,j)%ejcov + cumulative_elchange(i,j),0.0_DP)
      end do
   end do
return
end subroutine crater_soften

