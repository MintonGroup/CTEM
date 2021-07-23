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
   use module_ejecta
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
   integer(I4B) :: i,j,k,l,m,xpi,ypi,n,ntot,ioerr,inc,xi,xf,yi,yf
   integer(I4B) :: maxhits = 1
   real(DP) :: lamleft,dburial,lambda,Kbar,diam,radius,Area,avgejc,sf
   real(DP),dimension(domain%pnum) :: dN,lambda_regolith,Kbar_regolith,lambda_bedrock,Kbar_bedrock
   real(DP),dimension(2)    :: rn  
   real(DP) :: superlen,rayfrac
   type(cratertype) :: crater
   real(DP),dimension(:,:),allocatable :: ejdistribution
   integer(I4B),dimension(:,:),allocatable :: ejisray
   

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

      Kbar_bedrock(i)  = 0.84_DP * (0.5_DP * nflux(1,i))**4 ! Intrinsic degradation function
      Kbar_regolith(i) = 0.84_DP * (0.5_DP * nflux(2,i))**4
       
      if (user%dosoftening) then
         Kbar_bedrock(i) = Kbar_bedrock(i) + user%soften_factor * (0.5_DP * nflux(1,i))**(user%soften_slope) 
         Kbar_regolith(i) = Kbar_regolith(i) + user%soften_factor * (0.5_DP * nflux(2,i))**(user%soften_slope)
      end if

   end do

   do j = 1,user%gridsize
      do i = 1,user%gridsize
         do n = 1,ntot
            dburial = EXFAC * 0.5_DP * nflux(1,n)
            if (surf(i,j)%ejcov > dburial) then
               lambda = lambda_regolith(n)
               Kbar = Kbar_regolith(n)
               diam = nflux(2,n)
            else
               lambda = lambda_bedrock(n)
               Kbar = Kbar_bedrock(n)
               diam = nflux(1,n)
            end if 
            if (diam > domain%smallest_crater) exit
            k = util_poisson(lambda)
            kdiff(i,j) = kdiff(i,j) + k * Kbar / Area
         end do
      end do
   end do
   kdiff(0,0) = kdiff(user%gridsize,user%gridsize)
   kdiff(user%gridsize + 1,user%gridsize + 1) = kdiff(1,1)
   kdiff(0,:) = kdiff(user%gridsize,:)
   kdiff(:,0) = kdiff(:,user%gridsize)
   kdiff(user%gridsize + 1,:) = kdiff(1,:)
   kdiff(:,user%gridsize + 1) = kdiff(:,1)

   ! Now add the superdomain contribution to crater degradation 
   if (user%dosoftening) then
      crater%ejdis = user%gridsize * user%pix * 0.5_DP
      crater%continuous = crater%ejdis / DISEJB
      radius = (crater%continuous / 2.348_DP)**(1.0_DP / 1.006_DP)
      crater%frad = radius
      crater%rad = radius
      inc = nint(min(crater%ejdis / user%pix, crater%frad * user%ejecta_truncation / user%pix)) + 1
      crater%xl = user%gridsize * user%pix * 0.5_DP
      crater%yl = user%gridsize * user%pix * 0.5_DP
      crater%xlpx = nint(crater%xl / user%pix)
      crater%ylpx = nint(crater%yl / user%pix)
      allocate(ejdistribution(-inc:inc,-inc:inc))
      allocate(ejisray(-inc:inc,-inc:inc))
      call ejecta_ray_pattern(user,surf,crater,inc,-inc,inc,-inc,inc,ejdistribution)
      do j = -inc,inc
         do i = -inc,inc
            if (ejdistribution(i,j) > 0.0001_DP) then
               ejisray(i,j) = 1
            else
               ejisray(i,j) = 0
            end if
         end do
      end do
      rayfrac = PI * inc**2 / (1.0_DP * count(ejisray == 1))
      deallocate(ejisray,ejdistribution)

      avgejc = sum(surf%ejcov) / user%gridsize**2
      do l = 1,domain%pnum
         if ((PI * (user%soften_size * nflux(1,l) * 0.5_DP)**2 <= (user%gridsize)**2 ).and.&
             (PI * (user%soften_size * nflux(2,l) * 0.5_DP)**2 <= (user%gridsize)**2 )) cycle

         dburial = EXFAC * 0.5_DP * nflux(1,n)
         if (avgejc > dburial) then
            radius = nflux(2,l) * 0.5_DP
         else
            radius = nflux(1,l) * 0.5_DP
         end if

         dN(l) = nflux(3,l) * user%interval * finterval
         Area = min(PI * (user%soften_size * radius)**2 , 4 * PI * user%trad**2)
         Area = max(Area - (user%pix * user%gridsize)**2,0.0_DP)
         if (Area == 0.0_DP) cycle
         superlen = 0.5_DP * (sqrt(Area) - user%pix * user%gridsize)

         lambda = dN(l) * Area
         if (lambda == 0.0_DP) cycle
         k = util_poisson(lambda)
         do m = 1, k
            
            ! generate random crater on the superdomain
            ! Set up size of crater and ejecta blanket
            crater%frad = radius
            crater%rad = radius
            crater%continuous = 2.348_DP * crater%frad**(1.006_DP) 
            crater%ejdis = DISEJB * crater%continuous
            inc = nint(min(crater%ejdis / user%pix, crater%frad * user%ejecta_truncation / user%pix)) + 1

            ! find the x and y position of the crater
            call random_number(rn)
            crater%xl = real(superlen * (2 * rn(1) - 1.0_DP), kind=SP)
            crater%yl = real(superlen * (2 * rn(2) - 1.0_DP), kind=SP)

            if (crater%xl > 0.0_SP) then
               crater%xl = crater%xl + user%pix * user%gridsize
            else
               crater%xl = crater%xl - user%pix * user%gridsize
            end if
            if (crater%yl > 0.0_SP) then
               crater%yl = crater%yl + user%pix * user%gridsize
            else
               crater%yl = crater%yl - user%pix * user%gridsize
            end if

            crater%xlpx = nint(crater%xl / user%pix)
            crater%ylpx = nint(crater%yl / user%pix)

            xi = 1 - crater%xlpx
            xf = user%gridsize - crater%xlpx
            yi = 1 - crater%ylpx
            yf = user%gridsize - crater%ylpx

            allocate(ejdistribution(xi:xf,yi:yf))
            allocate(ejisray(xi:xf,yi:yf))
            ! Now generate ray pattern
            call ejecta_ray_pattern(user,surf,crater,inc,xi,xf,yi,yf,ejdistribution)
            do j = yi,yf
               do i = xi,xf
                  if (ejdistribution(i,j) > 0.0001_DP) then
                     ejisray(i,j) = 1
                  else
                     ejisray(i,j) = 0
                  end if
               end do
            end do

            do j = 1, user%gridsize
               do i = 1, user%gridsize
                  xpi = i - crater%xlpx 
                  ypi = j - crater%ylpx 
                  if ((abs(xpi) > inc) .or. (abs(ypi) > inc)) cycle
                  if (ejisray(xpi,ypi) == 0) cycle 
                  kdiff(i,j) = kdiff(i,j) + &
                         rayfrac * user%soften_factor / (PI * user%soften_size**2) * crater%frad**(user%soften_slope - 2.0_DP)
               end do
            end do
            deallocate(ejdistribution,ejisray)
         end do
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

