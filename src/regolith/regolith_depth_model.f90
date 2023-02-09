!**********************************************************************************************************************************
!
!  Unit Name   : regolith_mix
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Simulate a lunar regolith reworking zone by vertical mixing of several layers in our push-pop system        
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : surf : Surface expression matrix
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_depth_model(user,domain,finterval,nflux,p)
   use module_globals
   use module_util
   use, intrinsic :: ieee_arithmetic
   use module_regolith, EXCEPT_THIS_ONE => regolith_depth_model
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: finterval  ! time elapsed ratio to the total time 
   real(DP),dimension(:,:),intent(in) :: nflux ! impact rate (number of craters per m^2 per year)
   real(DP),dimension(:,:),intent(out) :: p

   ! Internal variables
   real(DP) :: a1, a2, rmax, rmin, dbin, psum, t
   real(DP),dimension(2,domain%pnum) :: nflux_pix
   integer(I4B) :: i, j
   real(DP)     :: h, dr, f, fmin, fmax
   real(DP)     :: ntotsubcrat
   logical :: underflow
   real(DP) :: psumfunc

   if (ieee_support_underflow_control(psum)) call ieee_set_underflow_mode(gradual=.true.)
   ! Smallest crater size in sub-pixel crater regime (regolith scaling column in "nflux")
   ! Do it in terms of number of craters
   ntotsubcrat = 0._DP
   do i = 1, domain%pnum
      ntotsubcrat = ntotsubcrat + nflux(3,i)
      nflux_pix(1,i) = nflux(1,i) / 2.0_DP / ( user%gridsize * user%pix )
      nflux_pix(2,i) = nflux(3,i) * domain%area * user%interval
   end do
   ntotsubcrat = ntotsubcrat * domain%area * user%interval
   nflux_pix(2,:) = nflux_pix(2,:) / ntotsubcrat
  
   t    = ntotsubcrat * finterval ! Time in unit of number of craters

   !do i = 1, domain%smallest_impactor_index
      !p(1,i) = nflux(1,i) * a 
      !h = nflux_pix(1,i) * 2.0 * a
      !psum = 0._DP
      !do j = i+1, domain%smallest_impactor_index-1
      !   rmax = nflux_pix(1,j)
      !   f    = nflux_pix(2,j) * (rmax**2 - h / (2.0 * a) * rmax)
      !   psum = psum + f 
      !end do
      !rmax = nflux_pix(1,domain%smallest_impactor_index)
      !f = 0.5 * nflux_pix(2,domain%smallest_impactor_index) * (rmax**2 - h / (2.0 * a) * rmax)
      !psum = psum + f
      !p(2,i) = 1.0_DP - exp(-1.0 * PI * psum * t)
   !end do
   !stop

   do i = 1, domain%smallest_impactor_index
      ! Loop over a wanted mixing depth
      p(1,i) = nflux(1,i) * ALPHA
      h      = nflux_pix(1,i) * 2.0_DP * ALPHA
      psum   = 0._DP
      do j = i, domain%smallest_impactor_index! Only craters 1/alpha times larger than mixing depth will be considered.
         ! Trapzoidal numerical integration of f(x) in [a, b]
         ! Non-uniform grid: sum up 0.5 * (r_i+1 - r_i) * (f(ri) + f(r_i+1)), where i in [a, b]
         ! f(x) here is dN(r)/dr * (r^2 - h/2a * r)
         ! f(ri) = dN(ri)/dri * (ri^2 - h/2a * ri)
         ! dN(ri) / dri, where dN(ri) is differential number of a crater's size, ri's bin in CTEM
         ! ri's bin in cumulative SFD (Nc) is | Nc(>r_i+1) - Nc(>ri) | * A * t with the bin 
         ! from ri to r_i+1. 
         ! So, dN(ri)/dri = (Nc(>ri) - Nc(>r_i+1)) * A * t / (r_i+1 - ri) 
         ! dN(r_i+1)/dr_i+1 = (Nc(>ri+1) -Nc(>r_i+2)) * A * t / (r_i+2 - ri) 
         ! So, f(r_i+1) can be calculated. 

         ! i-th bin of crater SFD
         rmax = nflux_pix(1,j+1)
         rmin = nflux_pix(1,j)
         dr   = rmax - rmin
         fmin = nflux_pix(2,j) / (nflux_pix(1,j+1) - nflux_pix(1,j)) * &
                ( rmin**2 - h / (2.0_DP * ALPHA) * rmin )
         fmax = nflux_pix(2,j+1) / (nflux_pix(1,j+2) - nflux_pix(1,j+1)) * &
                ( rmax**2 - h / (2.0_DP * ALPHA) * rmax )
         f    = 0.5_DP * dr * ( fmin + fmax )
         psum = psum + f   
      end do
      psumfunc = exp(-1.0_DP * PI * psum * t)
      if (psumfunc > epsilon(psum)) then
         p(2,i) = 1.0_DP - exp(-1.0_DP * PI * psum * t)
      else
         p(2,i) = 1.0_DP
      end if
   end do

   return
end subroutine regolith_depth_model
