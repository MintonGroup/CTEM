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
   use module_regolith, EXCEPT_THIS_ONE => regolith_depth_model
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: finterval  ! time elapsed ratio to the total time 
   real(DP),dimension(:,:),intent(in) :: nflux ! impact rate (number of craters per m^2 per year)
   real(DP),dimension(:,:),intent(out) :: p

   ! Internal variables
   real(DP), parameter :: a = 0.125
   real(DP) :: a1, a2, rmax, rmin, dbin, psum, t
   real(DP),dimension(2,domain%pnum) :: nflux_pix
   integer(I4B) :: i, j
   real(DP)     :: h, dr
   real(DP)     :: ntotsubcrat

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
  
   rmin = nflux_pix(1,1) 
   t    = ntotsubcrat * finterval ! Time in unit of number of craters

   do i = 1, domain%smallest_impactor_index
      p(1,i) = nflux(1,i) / 2.0_DP
      h = nflux_pix(1,i)
      psum = 0._DP
      do j = 1, domain%smallest_impactor_index
         if (nflux_pix(1,j) >= h / (2.0_DP * a)) then
            if (j==domain%smallest_impactor_index) then
               dr = ( ( log( nflux_pix(1,j) - nflux_pix(1,j-1) ) )**2 ) / &
                    ( log( nflux_pix(1,j-1) - nflux_pix(1,j-2) ) )
               rmax = nflux_pix(1,j) + exp(dr)            
               rmin = nflux_pix(1,j)
               dbin = exp(dr)
            else
               rmax = nflux_pix(1,j+1) 
               rmin = nflux_pix(1,j)
               dbin = rmax - rmin
            end if
            a1   = 1.0_DP/3.0_DP * rmax**3 - h / (4.0_DP * a) * rmax**2
            a2   = 1.0_DP/3.0_DP * rmin**3 - h / (4.0_DP * a) * rmin**2
            psum = psum + (-1.0_DP * PI * (a1 - a2) * nflux_pix(2,j) / dbin)
         end if
      end do

      p(2,i) = 1.0_DP - exp(psum * t)

   end do
   return
end subroutine regolith_depth_model
