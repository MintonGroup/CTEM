!****f* ejecta/ejecta_ray_pattern
! Name
!   ejecta_ray_pattern -- Calculate ejecta ray pattern
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_io
!   * module_crater
!   * module_regolith
!   * module_ejecta
!   
!   call ejecta_ray_pattern(user,surf,crater,inc,ejdistribution)
!
! DESCRIPTION
!    
!   Models the discontinuous ray pattern based on the Superformula. 
!   Citation: Gielis, J. "A Generic Geometric Transformation that Unifies a Wide Range of Natural 
!   and Abstract Shapes." Amer. J. Botany 90, 333-338, 2003. 
! 
! ARGUMENTS
!   Input
!   * user    -- User input parameters
!   * surf    -- Surface ggrid
!   * crater  -- Crater dimension container
!   * domain  -- Simulation domain variable container
!
!   Output
! 
!***

!**********************************************************************************************************************************
subroutine ejecta_ray_pattern(user,crater,i,j,diffi,eji)
   use module_globals
   use module_util
   use module_io
   use module_crater
   use module_regolith
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_ray_pattern
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   integer(I4B),intent(in) :: i,j
   real(DP),intent(out) :: diffi
   real(DP),intent(out) :: eji

   ! Internal variables
   integer(I4B) :: nrays,k,n,nef,incsq,iradsq,xpi,ypi,ejpxsq
   real(DP) :: frac,mef,lrad,lradsq,xp,yp,binres,areafrac,xbar,ybar
   real(DP) :: rn
   real(DP) :: theta, lradp, maxdistance
   real(DP), parameter :: n1 = 4.0_DP
   real(DP) :: n2, mag
   real(DP),dimension(xi:xf,yi:yf) :: isray
   real(DP),dimension(:),allocatable :: numinray,totnum
   real(DP),dimension(:),allocatable :: mefarray
   real(DP) :: ans
   logical :: bigej

   real(DP) :: C1,C2,p ! Fits to fe vs fd equation for the new ray pattern
   real(DP) :: rmin,rmax,r

   ! Flowery ray variables
   integer(I4B)  :: nfrays
   real(DP)      :: n1f, n2f, magf, lradf, rayf,rayavg
   real(DP), parameter :: a = 16.8799 !a = 11.8126 ! Fitting parameters for a relation between ray length and radius of crater 
   real(DP), parameter :: b = 0.120621 !0.143   ! based on Jake's crater rays mapping studies! 
   real(DP) :: mvrld                  ! median value of ray length distribution
   real(DP) :: mvrldsc                ! median value of ray length distribution scaled by continuous ejecta ext/allent
   logical :: ej
   real(DP),dimension(Nraymax) :: thetari

   eji = 0.0_DP
   diffi = 0.0_DP

   if (user%dorays) then
      do i = 1,Nraymax
         thetari(i) = 2 * pi * i / Nraymax
      end do
      call shuffle(thetari) ! randomize the ray pattern

      call random_number(rn) ! randomize the ray orientation
      rmax = user%ejecta_truncation 
      rmin = crater%continuous / crater%frad
      crater%fe = 10.0_DP ! Estimate the equivalent degradation radius
      xpi = crater%xlpx + i
      ypi = crater%ylpx + j

      ! Find distance from crater center to current pixel center in real space
      xp = xpi * user%pix
      yp = ypi * user%pix

      xbar = xp - crater%xl 
      ybar = yp - crater%yl

      areafrac = util_area_intersection(user%ejecta_truncation * crater%frad,xbar,ybar,user%pix) 
      r = sqrt(xbar**2 + ybar**2) / crater%frad
      theta = mod(atan2(ybar,xbar) + pi + rn * 2 * pi,2 * pi)
      diffi = areafrac * ejecta_ray_pattern_func(theta,r,rmin,rmax,thetari,.false.) 
      eji = areafrac * ejecta_ray_pattern_func(theta,r,rmin,rmax,thetari,.true.) 

   else
      !Do simple circular region
      incsq = inc**2
      iradsq = i*i + j*j

      if (iradsq < incsq) then

         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! Find distance from crater center to current pixel center in real space
         xp = xpi * user%pix
         yp = ypi * user%pix

         xbar = xp - crater%xl 
         ybar = yp - crater%yl
         areafrac = util_area_intersection(user%ejecta_truncation * crater%frad,xbar,ybar,user%pix) ! uniform circular
         diffi = areafrac
         eji = areafrac
      end if

   end if

   return  
   contains

   ! pure function ray(theta,thetar,r,n,w) result(ans)
   ! implicit none
   ! real(DP) :: ans
   ! real(DP),intent(in) :: theta,thetar,r,w
   ! integer(I4B),intent(in) :: n
   ! real(DP) :: thetap,thetapp,a,b,c,dtheta

   ! c = w / r
   ! b = thetar 
   ! dtheta = min(2*pi - abs(theta - b),abs(theta - b))
   ! a = sqrt(2 * pi) / (n * c * erf(pi / (2 *sqrt(2._DP) * c)))
   ! ans = a * exp(-dtheta**2 / (2 * c**2))

   ! return
   ! end function ray

   ! !pure function pattern(theta,r,rmin,rmax,thetari,ej) result(ans)
   ! !function pattern(theta,r,rmin,rmax,thetari,ej) result(ans)
   ! subroutine pattern(theta,r,rmin,rmax,thetari,ej,ans)
   ! implicit none
   ! real(DP) :: ans
   ! real(DP),intent(in) :: r,rmin,rmax,theta
   ! real(DP),dimension(:),intent(in) :: thetari
   ! logical,intent(in) :: ej
   ! real(DP) :: a,c
   ! real(DP) :: thetar,rw,rw0,rw1
   ! real(DP) :: f,rtrans,length,rpeak,minray,FF
   ! integer(I4B) :: n,i


   ! minray = rmin * 3

   ! if (r > rmax) then
   !    ans = 0._DP
   ! else if (r < 1.0_DP) then
   !    if (ej) then
   !       ans = 1.0_DP
   !    else
   !       ans = 0.0_DP
   !    end if
   ! else
   !    rw0 = rmin * pi / Nraymax / 2
   !    rw1 = 2 * pi / Nraymax
   !    rw = rw0 * (1._DP - (1.0_DP - rw1 / rw0) * exp(1._DP - (r / rmin)**2))
   !    n = max(min(floor((Nraymax**rayp - (Nraymax**rayp - 1) * log(r/minray) / log(rray/minray))**(1._DP/rayp)),Nraymax),1) ! Exponential decay of ray number with distance
   !    ans = 0._DP
   !    rtrans = r - 1.0_DP
   !    c = rw / r
   !    a = sqrt(2 * pi) / (n * c * erf(pi / (2 *sqrt(2._DP) * c)))
   !    do i = 1,Nraymax
   !       length = minray * exp(log(rray/minray) * ((Nraymax - i + 1)**rayp - 1_DP) / ((Nraymax**rayp - 1)))
   !       rpeak = (length - 1_DP) * 0.5_DP
   !       if (ej) then
   !          FF = 1.0_DP
   !          if (r > length) then
   !             f = 0.0_DP
   !          else
   !             f = a 
   !          end if
   !       else
   !          FF = rayfmult * (20 / rmax)**(0.5_DP) * 0.25_DP 
   !          f = FF * fpeak * (rtrans / rpeak)**rayq * exp(1._DP / rayq * (1.0_DP - (rtrans / rpeak)**rayq)) 
   !       end if
   !       ans = ans + ray(theta,thetari(i),r,n,rw) * f / a 
   !    end do
   ! end if 


   ! !end function pattern
   ! end subroutine pattern
   
   subroutine shuffle(a)
   real(DP), intent(inout) :: a(:)
   integer :: i, randpos
   real(DP) :: r,temp
    
   do i = size(a), 2, -1
     call random_number(r)
     randpos = int(r * i) + 1
     temp = a(randpos)
     a(randpos) = a(i)
     a(i) = temp
   end do

  end subroutine shuffle

    !^^^^^^^^^^^^^^^^^^^^^^^^

   ! *************************** Superformula Ray Model              ************************************!
   ! *************************** Part I.  Spoke and Skinny Ray Model ************************************!
   ! From fitting Jake Elliot's ray mapping data, it is a linear function that describes the relationship between the median value of ray length 
   ! distribution and the radius of rayed craters (in unit of kilometers)
   ! Also, we need to scale it with the continuous ejecta's extent for ray model
!#   mvrld      = a * (crater%frad/1000.0)**(b)
!#   mvrldsc    = mvrld / (continuous/crater%frad)
   ! It appears that no strong correlation between the number of rays and size of craters.
   ! The average number of rays is about 10. The minimum and maximum number is 6 and 14 respectively.

   ! Determine parameters of Ray model based on Gielis's superformula
!   ! There are three parameters regarding to our desired ray shape: n1, n2, m
   ! m:  the repeating part of formula, and it controls the number of rays (arms).
   ! n1: the default value is set as 4.0 in our rays case
!   ! n2: n2 and n1 combining together is to control the slenderness of a ray, in general, n1/n2 is always smaller than 1 in our cases.
!!   !     The smaller the ratio, the skinnier a ray.

!   nrays      = nint(8 * rn(1)) + 6
!   mvrld      = crater%ejdis
!   mvrldsc    = mvrld / crater%continuous
!   n2         = 8.0_DP * ( log10(mvrldsc/1.5_DP) / log10(2.0_DP) ) + 2.0_DP
!   ! *************************** Part II. Flowery Ray Model         ************************************!
!   nfrays = 20
!   n1f  = 1.0_DP
!   n2f  = 0.5_DP
!   rayf = 4.0 !7.5_DP
!
!   isray = 0.0_DP
!   incsq = inc**2
!   maxdistance = inc * user%pix

   ! Mass enhancement factor array
!   ! Don't let the bin resolution get too small for small craters. 
!   binres = min(10.0_DP,crater%ejdis / 10) * user%pix
!   nef = ceiling(crater%ejdis / binres) ! bin size for enhancement factor calculation
!   allocate(numinray(nef),totnum(nef),mefarray(nef))

!   totnum = 0
!   numinray = 0
!   do j = yi,yf
!      do i = xi,xf
!         iradsq = i*i + j*j

!         xpi = crater%xlpx + i
!         ypi = crater%ylpx + j

!         ! Find distance from crater center to current pixel center in real space
!         xp = xpi * user%pix
!         yp = ypi * user%pix
!         
!         lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2
!         lrad = sqrt(lradsq) 

!         if ((iradsq > incsq).or.(lrad <= crater%rad)) cycle
!
!         ! Sum up the number of pixels in this bin for the mass enhancement calculation
!         k = ceiling(lrad / binres)
!         totnum(k) = totnum(k) + 1.0_DP

!         ! Ray pattern setup
!         ! Take average of five points in the pixel to "feather" the edges of the ray
!         rayavg = 0.0_DP
!         do n = 1, 5
!            select case(n)
!            case(1)
!               theta = atan2(j * 1._DP,i * 1._DP) + 2.0_DP * PI
!            case(2)
!               theta = atan2(j * 1._DP - 0.5_DP,i * 1._DP - 0.5_DP) + 2.0_DP * PI
!               lradsq = (crater%xl - xp - 0.5_DP * user%pix)**2 + (crater%yl - yp - 0.5_DP * user%pix)**2
!               lrad = sqrt(lradsq) 
!            case(3)
!               theta = atan2(j * 1._DP - 0.5_DP,i * 1._DP + 0.5_DP) + 2.0_DP * PI
!               lradsq = (crater%xl - xp - 0.5_DP * user%pix)**2 + (crater%yl - yp + 0.5_DP * user%pix)**2
!               lrad = sqrt(lradsq) 
!            case(4)
!               theta = atan2(j * 1._DP + 0.5_DP,i * 1._DP + 0.5_DP) + 2.0_DP * PI
!               lradsq = (crater%xl - xp + 0.5_DP * user%pix)**2 + (crater%yl - yp + 0.5_DP * user%pix)**2
!               lrad = sqrt(lradsq) 
!            case(5)
!               theta = atan2(j * 1._DP + 0.5_DP,i * 1._DP - 0.5_DP) + 2.0_DP * PI
!               lradsq = (crater%xl - xp + 0.5_DP * user%pix)**2 + (crater%yl - yp - 0.5_DP * user%pix)**2
!               lrad = sqrt(lradsq) 
!            end select
!
!            mag   = ( ( (abs(cos(nrays * theta / 4.0_DP)))**n2 + &
!                    (abs(sin(nrays * theta / 4.0_DP)))**n2 )**(-1.0_DP/n1) ) 
!            magf  = rayf * ( ( (abs(cos(nfrays * theta / 4.0_DP)))**n2f + &
!                      (abs(sin(nfrays * theta / 4.0_DP)))**n2f )**(-1.0_DP/n1f))
!            lradp = crater%frad * mag
!            lradf = crater%frad * magf
!            lradp = lradp + lradf !max(lradp, lradf) 
!
!            if ((lrad < maxdistance) .and. (lrad < lradp)) then 
!               rayavg = rayavg + 1.0_DP
!            end if
!         end do
!         rayavg = rayavg / 5
!         isray(i,j) = rayavg
!         numinray(k) = numinray(k) + rayavg
!      end do
!   end do
!   if ((xi /= -inc).or.(xf /= inc).or.(yi /= -inc).or.(yf /= inc)) then
!      ejdistribution(xi:xf,yi:yf) = isray(xi:xf,yi:yf) 
!   else
!      do k = 1,nef
!         if (numinray(k) == 0) then
!            mefarray(k) = 1.0_DP
!         else
!            mefarray(k) = totnum(k) /  numinray(k)
!         end if
!      end do
!      ejdistribution = 1.0_DP

      ! Now calculate the ejecta distribution with mass enhancement factor
!      do j = yi,yf
!         do i = xi,xf
!            ejdistribution(i,j) = 0.0_DP
!            iradsq = i*i + j*j
!            lrad = sqrt(lradsq) 
!            if (lrad < crater%continuous) then
!               ejdistribution(i,j) = 1.0_DP
!            else
!               frac = lrad / binres
!               k = ceiling(frac)
!               ! Interpolate the mass enhancement factor between array points
!               if (k > 1) then
!                  mef = mefarray(k) + (frac - k * 1.0_DP) * (mefarray(k) - mefarray(k - 1))
!               else
!                  mef = mefarray(k) + (frac - k * 1.0_DP) * (mefarray(k + 1) - mefarray(k))
!               end if
!
   !            ejdistribution(i,j) = mef * isray(i,j)
   !         end if
   !       end do
   !    end do
   !end if

   !deallocate(numinray,totnum,mefarray) 

end subroutine ejecta_ray_pattern

