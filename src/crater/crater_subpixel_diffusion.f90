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
subroutine crater_subpixel_diffusion(user,surf,prod,nflux,domain,finterval,kdiffin)
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
   real(DP),dimension(:,:),intent(inout) :: kdiffin

   ! Internal variables
   real(DP),dimension(0:user%gridsize + 1,0:user%gridsize + 1) :: cumulative_elchange,kdiff
   integer(I4B),dimension(2,0:user%gridsize + 1,0:user%gridsize + 1) :: indarray
   integer(I4B) :: i,j,k,xpi,ypi,n,ntot
   integer(I4B) :: maxhits = 1
   real(DP) :: lamleft,dburial,lambda,kappat,diam,radius,Area,avgejc
   real(DP),dimension(domain%pnum) :: dN,lambda_regolith,kappat_regolith,lambda_bedrock,kappat_bedrock


   ! Create box for soften calculation (will be no bigger than the grid itself)
   do j = 0,user%gridsize + 1
      do i = 0,user%gridsize + 1
         xpi = i
         ypi = j
         call util_periodic(xpi,ypi,user%gridsize)
         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi
         kdiff(i,j) = kdiffin(xpi,ypi)
      end do
   end do

   ntot = 1
   ! calculate the subpixel diffusion probability function
   Area = user%pix**2

   do i = 1,domain%pnum
      if ((nflux(1,i) > domain%smallest_crater).and.(nflux(2,i) > domain%smallest_crater)) exit
      ntot = i
      dN(i) = nflux(3,i) * user%interval * finterval

      lambda_bedrock(i) = dN(i) * Area
      lambda_regolith(i) = dN(i) * Area
       
      if (user%dosoftening) then
         kappat_bedrock(i) = user%soften_factor * (0.5_DP * nflux(1,i))**(user%soften_slope)
         kappat_regolith(i) = user%soften_factor * (0.5_DP * nflux(2,i))**(user%soften_slope)
      end if

   end do

   do j = 0,user%gridsize + 1
      do i = 0,user%gridsize + 1
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
            kdiff(i,j) = kdiff(i,j) + k * kappat / Area
         end do
      end do
   end do

   ! Now add the superdomain contribution to crater degradation 
   if (user%dosoftening) then
      avgejc = sum(surf%ejcov) / user%gridsize**2
      do i = 1,domain%pnum
         if ((user%soften_size * nflux(1,i) <= user%gridsize / SQRT2 ).and.&
             (user%soften_size * nflux(2,i) <= user%gridsize / SQRT2 )) cycle

         dburial = EXFAC * 0.5_DP * nflux(1,n)
         if (avgejc > dburial) then
            radius = nflux(2,i) * 0.5_DP
         else
            radius = nflux(1,i) * 0.5_DP
         end if

         dN(i) = nflux(3,i) * user%interval * finterval
         Area = min(2 * (user%soften_size * radius)**2 , PI * user%trad**2)
         Area = max(Area - user%gridsize**2,0.0_DP)

         lambda = dN(i) * Area
         if (lambda == 0.0_DP) cycle
         k = util_poisson(lambda)
          
         if (k > 0) kdiff = kdiff + k * user%soften_factor / (PI * user%soften_size**2) * radius**(user%soften_slope - 2.0_DP)
      end do
   end if
   

   call util_diffusion_solver(user,surf,user%gridsize + 2,indarray,kdiff,cumulative_elchange,maxhits)
   do j = 1,user%gridsize
      do i = 1,user%gridsize
         surf(i,j)%dem = surf(i,j)%dem + cumulative_elchange(i,j)
         surf(i,j)%ejcov = max(surf(i,j)%ejcov + cumulative_elchange(i,j),0.0_DP)
      end do
   end do
   kdiffin = 0.0_DP

return
end subroutine crater_subpixel_diffusion

