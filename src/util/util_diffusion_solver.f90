!**********************************************************************************************************************************
!
!  Unit Name   : util_diffusion_solver
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Solves the diffusion equation with spatially varying diffusion coefficient
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       : Uses a 2nd order central difference explicit finite difference scheme
!                Dirichlet boundary conditions assumed, such that cumulative_elchange at boundaries is 0
!
!**********************************************************************************************************************************
subroutine util_diffusion_solver(user,surf,N,indarray,kdiff,cumulative_elchange,maxhits)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_diffusion_solver
   use, intrinsic :: ieee_exceptions
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   integer(I4B),intent(in) :: N
   real(DP),dimension(N,N),intent(in) :: kdiff
   integer(I4B),dimension(2,N,N),intent(in) :: indarray
   real(DP),dimension(N,N),intent(out) :: cumulative_elchange 
   integer(I4B),intent(in) :: maxhits

   ! Internal variables
   integer(I4B) :: i,j,loop
   real(DP) :: kdiffmax,dt,tmax,fac
   real(DP) :: grad,xijp1,xijm1,xip1j,xim1j,xij,kij,kip1j,kim1j,kijp1,kijm1
   integer(I4B) :: mi,mip1,mim1,mj,mjp1,mjm1,nloops
   real(DP),dimension(N,N) :: elchange
   integer(I4B),dimension(6,N,N) :: indarray_extended
   logical, dimension(size(IEEE_ALL))   :: fpe_halting_modes, fpe_quiet_modes
   logical, dimension(size(IEEE_USUAL)) :: fpe_flag 

   call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
   fpe_quiet_modes(:) = .false.
   call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

   cumulative_elchange = 0.0_DP
   kdiffmax = maxval(kdiff)

   if (abs(kdiffmax) < VSMALL) return
   dt = 0.25_DP * user%pix**2 / kdiffmax / real(maxhits,kind=DP) ! Stability limit

   tmax = 1.0_DP ! Equations are assumed to be derived for a unit system where t = 1
   nloops = ceiling(tmax / dt)
   dt = tmax / nloops
   ! Make extended index array with points +1 and -1 away in x and y from the origin point at i,j
   ! This is done to reduce the expense of calling util_periodic every loop
   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(2*N > INCPAR) &
   !$OMP SHARED(user,N,indarray,indarray_extended) 
   do j = 1,N 
      do i = 1,N
         mi = indarray(1,i,j)
         mj = indarray(2,i,j)
         mip1 = mi + 1
         mjp1 = mj + 1
         mim1 = mi - 1
         mjm1 = mj - 1

         ! periodic boundary conditions 
         call util_periodic(mi,mj,user%gridsize)
         call util_periodic(mip1,mjp1,user%gridsize)
         call util_periodic(mim1,mjm1,user%gridsize)
         indarray_extended(1,i,j) = mi
         indarray_extended(2,i,j) = mj
         indarray_extended(3,i,j) = mip1
         indarray_extended(4,i,j) = mjp1
         indarray_extended(5,i,j) = mim1
         indarray_extended(6,i,j) = mjm1
      end do
   end do
   !$OMP END PARALLEL DO


   fac = 0.5_DP * dt / user%pix**2 ! Constant factor multiplied by every element on the grid
   !  begin soften diffusion loop(s)
   do loop = 1,nloops
      elchange = 0.0_DP
      ! loop over affected matrix area
      !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(2*N > INCPAR) &
      !$OMP SHARED(surf,N,elchange,cumulative_elchange,indarray_extended,kdiff,fac) 
      do j = 2, N - 1 
         do i = 2,N - 1 
            ! Retrieve coordinates and diffusion constants
            mi   = indarray_extended(1,i,j)
            mj   = indarray_extended(2,i,j)
            mip1 = indarray_extended(3,i,j)
            mjp1 = indarray_extended(4,i,j)
            mim1 = indarray_extended(5,i,j)
            mjm1 = indarray_extended(6,i,j)

            kij   = kdiff(i  ,j  ) 
            kip1j = kdiff(i+1,j  )
            kim1j = kdiff(i-1,j  ) 
            kijp1 = kdiff(i  ,j+1)
            kijm1 = kdiff(i  ,j-1)

            ! Second order central difference with spatially varying diffusivity
            xij   = surf(mi  ,mj  )%dem + cumulative_elchange(i  ,j  )
            xip1j = surf(mip1,mj  )%dem + cumulative_elchange(i+1,j  )
            xim1j = surf(mim1,mj  )%dem + cumulative_elchange(i-1,j  )
            xijp1 = surf(mi  ,mjp1)%dem + cumulative_elchange(i  ,j+1)
            xijm1 = surf(mi  ,mjm1)%dem + cumulative_elchange(i  ,j-1)

            grad  =        ( kij + kip1j ) * ( xip1j - xij   )  
            grad  = grad - ( kij + kim1j ) * ( xij   - xim1j )
            grad  = grad + ( kij + kijp1 ) * ( xijp1 - xij   )
            grad  = grad - ( kij + kijm1 ) * ( xij   - xijm1 )

            elchange(i,j) = grad * fac
         end do
      end do 
      !$OMP END PARALLEL DO
      ! add the total diffusion to the cumulative
      cumulative_elchange = cumulative_elchange + elchange 
   end do

   call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Restore the original FPE halting modes

return
end subroutine util_diffusion_solver


