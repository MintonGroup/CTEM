!**********************************************************************************************************************************
!
!  Unit Name   : init_domain
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the simulation domain
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine init_domain(user,crater,domain,prod,pdist,vdist,crtscl)
   use module_globals
   use module_crater
   use module_ejecta
   use module_init, EXCEPT_THIS_ONE => init_domain
   implicit none

   ! Arguments
   type(usertype),intent(in)                       :: user
   type(cratertype),intent(inout)                  :: crater
   type(domaintype),intent(inout)                  :: domain
   real(DP),dimension(:,:),intent(inout)           :: prod,vdist
   real(DP),dimension(:,:),intent(out)             :: crtscl
   real(DP),dimension(:,:),intent(out),allocatable :: pdist

   ! Internals
   integer(I4B),parameter :: tistfac=100 ! Temporary storage array for crater tally is 100x larger than the displayed bin size
   real(DP) :: fcrat,bedrockfcrat,regolithfcrat,rmsvel,disthi,diffnum
   integer(I4B) :: i,j,k,plo,p
   real(DP) :: bedrockejdis,regolithejdis,ejdis
   real(DP) :: numvel,vlo
   type(ejbtype),dimension(EJBTABSIZE) :: ejb_r,ejb_b,ejb
   integer(I4B) :: ejtble_r,ejtble_b,ejtble,inc
   real(DP) :: h,hmean,hmeanbig,area,r1,r2,f1,f2,volume,hprofile,rmax,lrad,lradsq,newelev,melev,vtot
   type(surftype),dimension(0:EJBTABSIZE) :: profilesurf
   logical :: firstrun
    
   ! Executable code

   ! Some preliminary calculations to set the size of the domain
   domain%side = user%pix * user%gridsize
   domain%GM = user%gaccel * user%trad**2
   domain%parea = user%pix**2
   domain%area = domain%side**2
   domain%biggest_crater = domain%side * user%maxcrat
   domain%smallest_crater = user%pix
   domain%smallest_ejecta = SMALLESTEJECTA * user%pix
   domain%smallest_counted_crater = 2._DP / (1._DP + COUNTINGRIM) * sqrt(SMALLESTCOUNTABLE / PI) * user%pix 
   domain%vescsq = 2 * user%gaccel * user%trad


   ! Set up transition values
   select case(user%mat)
   case("ROCK")
      crater%cxexp = CXEXPS
      crater%cxtran = SIMCOMKS*(user%gaccel**SIMCOMPS)
   case("ICE")
      crater%cxexp = CXEXPI
      crater%cxtran = SIMCOMKI*(user%gaccel**SIMCOMPI)
   end select

   ! Now we build an idealized production population
   domain%initialize = .true.
   if (user%doangle) then 
      crater%sinimpang = 0.5_DP * SQRT2 ! Use 45 degree impact angle if user is allowing impact angle to vary
   else
      crater%sinimpang = 1.0_DP ! Use 90 degree impact angle if user is not allowing impact angle to vary
   end if

   rmsvel = 0._DP
   numvel = 0._DP
   vlo = -1.0_DP
   do k=1,domain%vnum
      if ((vlo < 0.0_DP).and.(vdist(2,k) > 0._DP)) then
         vlo = vdist(1,k)
         domain%vlo = k
      end if
      if (vdist(3,k) < 1.0_DP) domain%vhi = k 
      rmsvel = rmsvel + vdist(2,k)*vdist(1,k)**2
      numvel = numvel + vdist(2,k)
   end do
   domain%vhi = min(domain%vhi + 1,domain%vnum)
   rmsvel = sqrt(rmsvel/numvel) 

   crater%strflag = 1 ! Set to bedrock value
   crater%impvel = rmsvel
   do k=1,domain%pnum
      crater%imp  = prod(1,k)
      call crater_generate(user,crater,domain)
      crtscl(1,k) = crater%imp
      crtscl(2,k) = crater%fcrat
      ! Build a logspace version of the table for more efficient calculations later
      prod(3,k) = log(prod(1,k))
      prod(4,k) = log(prod(2,k))
   end do

   ! Bin the production distribution
   plo=100 ! To be safe, we begin at a crater that is 7500 AU wide
   do 
      disthi=1e3_DP*SQRT2**plo
      if (disthi > crtscl(2,domain%pnum)) then
         domain%pdistl = plo + 1
      end if 
      if (disthi < crtscl(2,1)) then
         exit
      else
         plo = plo - 1
      end if
   end do
   domain%pdistl = domain%pdistl - plo

   allocate(pdist(6,domain%pdistl))
   do i = 1,domain%pdistl
      pdist(1,i) = 1e3_DP*SQRT2**(plo+(i-1))
      pdist(2,i) = 1e3_DP*SQRT2**(plo+i)
      pdist(3,i) = 0.0_DP 
      pdist(4,i) = 0.0_DP
      pdist(5,i) = 0.0_DP
      pdist(6,i) = 0.0_DP
   end do

   do k = 1,domain%pnum
      i = ceiling(log(crtscl(2,k)/1e3_DP)/LOGSQRT2) - plo
      if (k == domain%pnum) then
         diffnum = prod(2,k)
      else
         diffnum = prod(2,k) - prod(2,k+1)
      end if
      pdist(3,i) = pdist(3,i) + diffnum * log(crtscl(2,k)) ! Geometric mean (intermediate step)
      pdist(4,i) = pdist(4,i) + diffnum      ! Differential number
   end do
   

   do i = 1,domain%pdistl
      if (pdist(4,i) > 0._DP) then
         pdist(3,i) = exp(pdist(3,i) / pdist(4,i)) ! Geometric mean (final step)
      else
         pdist(3,i) = sqrt(pdist(1,i) * pdist(2,i))
      end if
      pdist(5,i) = sum(pdist(4,i:domain%pdistl)) ! Cumulative number
      pdist(6,i) = (pdist(4,i)*pdist(3,i)**3)/(domain%area*(pdist(2,i)-pdist(1,i))) ! R-value
   end do
   

   ! Next we determine the smallest impactor that we ever expect to consider in this simulation
   ! This will be the smallest of: the impactor that produces a crater at least 1 pixel wide at the maximum velocity, the impactor
   ! that produces an ejecta blanket 3 pixels wide, or the smallest impactor in the production SFD
   
   domain%smallest_impactor_index = 1
   ! Find the smallest impactor that produces a crater at least 1 pixel wide at the maximum possible impact velocity and angle
   crater%impvel = vdist(1,domain%vnum)
   crater%sinimpang = 1.0_DP
   do k=1,domain%pnum
      crater%imp  = prod(1,k)

      ! Calculate crater size based on both strength models
      crater%strflag = 0
      call crater_generate(user,crater,domain)
      regolithfcrat = crater%fcrat

      crater%strflag = 1
      call crater_generate(user,crater,domain)
      bedrockfcrat = crater%fcrat

      fcrat = max(regolithfcrat,bedrockfcrat)
      domain%subcrater_limit = min(regolithfcrat,bedrockfcrat)
      if (fcrat >= user%pix) then
         domain%smallest_impactor_index = k
         exit
      endif
   end do

   goto 1000 
   ! Now estimate the amount of topographic overturn produced by subpixel craters
   ! This will be modeled as a diffusion rate
   crater%impvel = rmsvel
   k = domain%smallest_impactor_index
   crater%sinimpang = 0.5_DP
   domain%subpixel_ejecta_thickness = 0.0_DP
   call ejecta_distance_estimate(user,crater,domain,ejdis)
   vtot = 0._DP
   hmean = 0._DP
   hmeanbig = 0._DP
   do k = 1,domain%pnum 
   !do k = domain%smallest_impactor_index,1,-1
      crater%imp  = prod(1,k)

      ! Calculate crater size based on both strength models
      crater%strflag = 0
      call crater_generate(user,crater,domain)
      call ejecta_table_define(user,crater,domain,ejb,ejtble)
      ejdis = crater%ejdis
      if (crater%fcrat > domain%biggest_crater) exit

      area = 0.0_DP
      volume = 0.0_DP

      ! Estimate the total volume of ejecta produced per impact
      do i = 2,ejtble 
         r1 = exp(ejb(i-1)%lrad)
         r2 = exp(ejb(i)%lrad)
         h = exp(ejb(i-1)%thick)
         volume = volume +  PI *  (r2**2 - r1**2) * h
      end do
     
      ! Now estimate the volume of material excavated by the craterform itself
      ! We'll construct a crater profile
      !domain%small = crater%frad * SMALLFAC
      !call crater_find_visible(user,crater,domain)
      !rmax = crater%frad * (domain%small/crater%rheight)**(-1._DP/RIMDROP) !  Maximum distance of crater form
      !profilesurf%dem = 0._DP
      !newelev = 0._DP
      !melev = 0._DP
      !do i = 0,EJBTABSIZE
      !   lrad = rmax * i / real(EJBTABSIZE,kind=DP)
      !   lradsq = lrad**2
      !   if (lrad < crater%frad) then
      !      call crater_form_interior(user,profilesurf(i),crater,lradsq,newelev,melev)
      !   else 
      !      call crater_form_exterior(profilesurf(i),crater,domain,lradsq,newelev) 
      !   end if 
      !end do
      !do i = 1,EJBTABSIZE
      !   r1 = rmax * (i - 1) / real(EJBTABSIZE,kind=DP)
      !   r2 = rmax * i  / real(EJBTABSIZE,kind=DP)
      !   f1 = (profilesurf(i - 1)%dem)
      !   f2 = (profilesurf(i)%dem)
      !   h = 0.5_DP * (f1 + f2)
      !   volume = volume +  PI *  (r2**2 - r1**2) * abs(h)
      !end do
      if (k < domain%smallest_impactor_index) then
         hmean = hmean + volume * (prod(2, k) - prod(2,k + 1)) 
      else
         hmeanbig = hmeanbig + volume * (prod(2, k) - prod(2,k + 1)) 
      end if
      !domain%subpixel_ejecta_thickness = max(domain%subpixel_ejecta_thickness, hmean )
   end do
   !write(*,*) "hmean / hmeanbig = ",hmean/hmeanbig
   !domain%subpixel_ejecta_thickness = hmean * 2 * user%pix / domain%area
   !write(*,*) "Smallest impacts per unit area: ",(prod(2,1) - prod(2,2)) / domain%area
   !write(*,*) "Total subpixel diffusion equivalent thickness: ",domain%subpixel_ejecta_thickness
   !domain%subpixel_ejecta_thickness = user%sf 
   !write(*,*) "User-supplied diffusion equivalent thickness:  ",domain%subpixel_ejecta_thickness
   !domain%subpixel_ejecta_thickness =  domain%subpixel_ejecta_thickness / prod(2,domain%smallest_impactor_index)
1000 continue
   ! Find the lowest value for the distribution (pinned at 1 km)
   p=100 ! To be safe, we begin at a crater that is 7500 AU wide
   do 
      disthi=1e3_DP*SQRT2**p 
      if (disthi < domain%subcrater_limit) then
         exit
      else
         p = p - 1
      end if
   end do
   !p = p + 1 ! Ignore the lowest bin because it may straddle the lower limit 
   domain%plo = p

   ! Find the largest value for the distribution
   domain%distl = 1
   do
      disthi=1e3_DP*SQRT2**p
      if (disthi > domain%side) exit
      domain%distl = domain%distl+1
      p = p + 1
   end do

   domain%initialize = .false.
   domain%small = user%pix * SMALLFAC
   return

end subroutine init_domain
