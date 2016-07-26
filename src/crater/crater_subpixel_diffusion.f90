!**********************************************************************************************************************************
!
!  Unit Name   : crater_subcrater_diffusion
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
subroutine crater_subpixel_diffusion(user,surf,prod,nflux,domain,finterval)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_subpixel_diffusion
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   real(DP),dimension(:,:),intent(in) :: prod,nflux 
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: finterval

   ! Internal variables
   real(DP),dimension(0:user%gridsize + 1,0:user%gridsize + 1) :: cumulative_elchange,kdiff
   integer(I4B),dimension(2,0:user%gridsize + 1,0:user%gridsize + 1) :: indarray
   integer(I4B) :: i,j,k,xpi,ypi,n,ntot
   integer(I4B) :: maxhits = 1
   real(DP) :: lamleft,dburial,lambda,kappat,diam
   real(DP),dimension(domain%pnum) :: dN,lambda_regolith,kappat_regolith,lambda_bedrock,kappat_bedrock

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
  
   ntot = 1
   ! calculate the subpixel diffusion probability function
   do i = 1,domain%pnum
      if ((nflux(1,i) > domain%smallest_crater).and.(nflux(2,i) > domain%smallest_crater)) exit
      ntot = i
      dN(i) = nflux(3,i) * user%interval * finterval

      lambda_bedrock(i) = dN(i) * 0.25_DP * PI * nflux(1,i)**2 
      lambda_regolith(i) = dN(i) * 0.25_DP * PI * nflux(2,i)**2 

      if (user%dosoftening) then
         kappat_bedrock(i) = SOFTEN_FACTOR * nflux(1,i)**(SOFTEN_SLOPE)
         kappat_regolith(i) = SOFTEN_FACTOR * nflux(2,i)**(SOFTEN_SLOPE)
      end if

   end do

   kdiff = 0._DP
   do j = 1,user%gridsize
      do i = 1,user%gridsize
         do n = 1,ntot
            dburial = EXFAC * 0.5_DP * nflux(1,n)
            if (surf(i,j)%ejcov > dburial) then
               lambda = lambda_regolith(n)
               kappat = kappat_regolith(n)
               diam = nflux(2,n)
            else
               lambda = lambda_bedrock(n)
               kappat = kappat_bedrock(n)
               diam = nflux(1,n)
            end if 
            if (diam > domain%smallest_crater) exit
            k = util_poisson(lambda)
            kdiff(i,j) = kdiff(i,j) + k * kappat
         end do
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
end subroutine crater_subpixel_diffusion

