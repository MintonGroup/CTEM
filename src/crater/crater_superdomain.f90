!**********************************************************************************************************************************
!
!  Unit Name   : crater_superdomain
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
subroutine crater_superdomain(user,surf,age,age_resolution,prod,nflux,domain,finterval)
   use module_globals
   use module_util
   use module_ejecta
   use module_regolith
   use module_crater, EXCEPT_THIS_ONE => crater_superdomain
   implicit none

   ! Arguments
   type(usertype),intent(in)                           :: user
   type(surftype),dimension(:,:),intent(inout)         :: surf
   real(DP),intent(in)                                 :: age
   real(DP),intent(in)                                 :: age_resolution
   real(DP),dimension(:,:),intent(in)                  :: prod,nflux 
   type(domaintype),intent(in)                         :: domain
   real(DP),intent(in)                                 :: finterval

   ! Internal variables
   integer(I4B),dimension(2,0:user%gridsize + 1,0:user%gridsize + 1) :: indarray
   integer(I4B) :: i,j,k,l,m,xpi,ypi,n,ntot,ioerr,inc,xi,xf,yi,yf
   real(DP) :: lamleft,dburial,lambda,diam,radius,Area,avgejc,sf
   real(DP),dimension(domain%pnum) :: dN
   real(DP),dimension(2)    :: rn  
   real(DP) :: superlen,rayfrac
   type(cratertype) :: crater
   real(DP),dimension(:,:),allocatable :: ejdistribution
   integer(I4B),dimension(:,:),allocatable :: ejisray

   ! Melt or glassy ray test
   real(DP) :: xp, yp, lradsq, lrad, erad
   real(DP),dimension(4) :: lradif
   real(DP) :: rm, depthb, maxgcrat, maxgcratkm
   real(SP) :: gglass, glrad

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

   ! calculate the smallest crater size for a crater forming in regolith
   do i = 1,domain%pnum
      if ((nflux(1,i) > domain%smallest_crater).and.(nflux(2,i) > domain%smallest_crater)) exit
      ntot = i
   end do
   n = ntot
   avgejc = sum(surf%ejcov) / user%gridsize**2
   maxgcrat = 0.0_DP
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
      ! Center at the local domain, and calculate the maximum distance from the super
      ! domain crater  to the local domain. 
      superlen = 0.5_DP * (sqrt(Area) - user%pix * user%gridsize)
      if (superlen < 0.0_DP) cycle

      lambda = dN(l) * Area
      if (lambda == 0.0_DP) cycle
      k = util_poisson(lambda)
      maxgcratkm = 0.0 
      do m = 1, k
            
         ! generate random crater on the superdomain
         ! Set up size of crater and ejecta blanket
         crater%frad = radius
         crater%rad  = radius
         crater%continuous = 2.348_DP * crater%frad**(1.006_DP) 
         crater%ejdis = DISEJB * crater%continuous
         inc = nint(min(crater%ejdis / user%pix, crater%frad * DISEJB / user%pix)) + 1

         ! find the x and y position of the crater
         ! When a superdomain crater's center is at the edge of a space that
         ! superlen defines, it has least ejecta delivering to the local domain.
         ! Soften_size should be the maximum ray extennt (beyond visible rays)
         ! to account for cold and hot ejecta. 
         ! Note for age data: new stream lines overlapped with melts needs to
         ! calculate here, and one can check out regolith_streamtube: lines 251 - 292.  
         call random_number(rn)
         crater%xl = real(superlen * (2 * rn(1) - 1.0_DP), kind=SP)
         crater%yl = real(superlen * (2 * rn(2) - 1.0_DP), kind=SP)

         ! Superdoamin crater can be from any of four quardants that is space defined by
         ! superlen.
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

         ! Define the pixel space in the local domain coordinate. 
         xi = 1 - crater%xlpx
         xf = user%gridsize - crater%xlpx
         yi = 1 - crater%ylpx
         yf = user%gridsize - crater%ylpx

         lradif(1) = sqrt( (xi * user%pix) **2 + (yi * user%pix)**2 )
         lradif(2) = sqrt( (xi * user%pix) **2 + (yf * user%pix)**2 )
         lradif(3) = sqrt( (xf * user%pix) **2 + (yi * user%pix)**2 )
         lradif(4) = sqrt( (xf * user%pix) **2 + (yf * user%pix)**2 )
         lrad  = minval(lradif)
         if (lrad > crater%ejdis) cycle

         allocate(ejdistribution(xi:xf,yi:yf))
         allocate(ejisray(xi:xf,yi:yf))
         ! Now generate ray pattern
         call ejecta_ray_pattern(user,surf,crater,inc,xi,xf,yi,yf,ejdistribution)
         ! Now if doregotrack is on, do melt zone calculation
         if (user%doregotrack) call regolith_melt_zone_superdomain(user,crater,domain,rm,depthb)

         do j = yi,yf
            do i = xi,xf
               ! Ejecta ray distribution
               if (ejdistribution(i,j) > 0.0001_DP) then
                  ejisray(i,j) = 1
               else
                  ejisray(i,j) = 0
               end if  
            end do
         end do

         if (user%doregotrack) then
            do j = 1, user%gridsize
               do i = 1, user%gridsize
                  xpi = i - crater%xlpx 
                  ypi = j - crater%ylpx 
                  if ((abs(xpi) > inc) .or. (abs(ypi) > inc)) cycle
                  if (ejisray(xpi,ypi) == 0) cycle  
                  call regolith_superdomain(user,crater,domain,surf(i,j)%regolayer,ejdistribution(xpi,ypi),&
                       i,j,age,age_resolution,rm,depthb)
               end do
            end do

         end if

         deallocate(ejdistribution,ejisray)

       end do

   end do

   return
end subroutine crater_superdomain

