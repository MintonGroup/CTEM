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
subroutine init_domain(user,crater,domain,prod,pdist,vdist,crtscl,nflux)
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
   real(DP),dimension(:,:),intent(inout),allocatable, optional :: nflux

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

   !real(DP) :: ejdis_regolith, ejdis_bedrock
   ! Test sub-pixel mixing:
   real(DP) :: f, a_crat, t !f: the fraction of a surface disturbed by craters
                                   !a_crat: the area of a crater, which is pi * r^2. 
                                   !t: time
   real(DP) :: kappaNdot,dN
    
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

   !# Next, we determine the impactor that we will only generate the ejecta rather than an actual crater based on a limit of vertical mixing.
   !# Given that Gault's study on turnover depth on the Moon, the upper 0.5 mm layer would have turnovered more than 100 times in the first billion
   !# years. And the upper 1 cm deep layer would turnover one time in the first 10^7 years. If we apply Gault's turnover model on our vertical mixing
   !# , and we won't be bothered by enormous and astronomical number of small craters slowing down CTEM, and simply just want to simulate the effect 
   !# without making any other contribution to, for example, regolith growth. So, I suggest that we can take this turnover depth, 1 cm in 10^7 years, as
   !# our limit of craters that contribute to surface, and the size of this crater in diameter is eight times this turnover depth: 8 cm.   
   !domain%smallest_ejecta_index = 1
   !crater%impvel = vdist(1,domain%vnum)
   !crater%sinimpang = 1.0_DP
   !do k=1,domain%pnum
   !   crater%imp = prod(1,k)
      ! Do it as the above code trying to figure out the size of a crater with both strength models
   !   crater%strflag = 0 
   !   call crater_generate(user,crater,domain)
   !   regolithfcrat = crater%fcrat
   !   ejdis_regolith = 25.0 * 2.3_DP * crater%frad**(1.006_DP)
      
   !   crater%strflag = 1
   !   call crater_generate(user,crater,domain)
   !   bedrockfcrat = crater%fcrat
   !   ejdis_bedrock = 25.0 * 2.3_DP * crater%frad**(1.006_DP)
 
   !   fcrat = max(regolithfcrat, bedrockfcrat)
   !   ejdis = max(ejdis_regolith, ejdis_bedrock)
   !   if (ejdis > domain%smallest_ejecta) then
   !      domain%smallest_ejecta_index = k
   !      domain%smallest_ejecta_crater = fcrat
         !write(*,*) domain%smallest_ejecta_index, prod(1,domain%smallest_ejecta_index),&
         !domain%smallest_ejecta_crater
   !      exit
   !   end if
   !end do
   !write(*,*) domain%smallest_ejecta_index, fcrat

   ! Testing retriving the impactor production SFD for all craters.
   !write(*,*) domain%pnum
   !do i = 1,domain%pnum !smallest_impactor_index
   !diffnum = prod(2,i) - prod(2,i+1)
   !write(*,*) prod(1,i), diffnum
   !write(*,*) pdist(1,i), pdist(1,4)
   !end do

!   goto 1000 
   ! Now estimate the amount of topographic overturn produced by subpixel craters
   ! This will be modeled as a diffusion rate

   if (present(nflux)) then 
   allocate(nflux(2,domain%smallest_impactor_index))
   crater%impvel = rmsvel
   crater%sinimpang = 0.5_DP * SQRT2
   crater%strflag = 0
   do k = 1,domain%smallest_impactor_index
      crater%imp  = prod(1,k)
      ! Calculate crater size based on both strength models
      diffnum = prod(2,k) - prod(2,k+1)
      call crater_generate(user,crater,domain)
      !write(*,*) k, crater%fcrat, diffnum
      ! My-T model for turnover 
      ! Since we obtain the total differential number for a given size of a crater, 
      ! it allows us to calculate the average impact rate for a given size of a crater,
      ! which is defined as variable, "nflux". 
      ! It is, nflux = dN(k,k+1) / A_tot / t_tot.
      ! Then, we can use this to calculate the fraction of a area for a given time by
      ! a tanh function in whcih Toshi derived an analytical model, which is that the 
      ! probabiluty of filling the surface by the same sized circle is past-dependent 
      ! event. That is, based on the limited space, as the disturbed surface becomes larger, 
      ! the probability of a new circle overlapping with the previous disturbed area is 
      ! dependent on the surface area of the  disturbed area. 
      ! This function is: p = (1 - exp(-2*nflux*pi*r^2)) / (1 + exp(-2*nflux*pi*r^2)),
      ! where the r is the radius of a crater. 
      nflux(1,k) = crater%fcrat
      nflux(2,k) = diffnum/(user%gridsize * user%pix)**2/user%interval
      !a_crat = PI * (crater%frad)**2
      !t    = user%interval 
      !f(k) = (1.0_DP - exp(-2.0_DP * nflux * a_crat * t)) / (1.0_DP + exp(-2.0_DP * nflux * a_crat * t))
      !write(*,*) 2.0 * crater%frad, nflux, f
   end do
   end if

   !write(*,*) "hmean / hmeanbig = ",hmean/hmeanbig
   !domain%subpixel_ejecta_thickness = hmean * 2 * user%pix / domain%area
   !write(*,*) "Smallest impacts per unit area: ",(prod(2,1) - prod(2,2)) / domain%area
   !write(*,*) "Total subpixel diffusion equivalent thickness: ",domain%subpixel_ejecta_thickness
   !domain%subpixel_ejecta_thickness = user%sf 
   !write(*,*) "User-supplied diffusion equivalent thickness:  ",domain%subpixel_ejecta_thickness
   !domain%subpixel_ejecta_thickness =  domain%subpixel_ejecta_thickness / prod(2,domain%smallest_impactor_index)
!1000 continue
   ! Find the lowest value for the distribution (pinned at 1 km)
   p=100 ! To be safe, we begin at a crater that is 7500 AU wide
   do 
      disthi=1e3_DP*SQRT2**p 
      !if (disthi < domain%subcrater_limit) then
      if (disthi < domain%smallest_ejecta_crater) then ! Change the size of a bin in true crataer distribution for
                                                       ! accommodating smallest craters that generate pixeled ejecta
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

   domain%small = user%pix * SMALLFAC


   !Calculate the sub-pixel topographic diffusion
   !domain%subpixel_diffusion_const = 0.0_DP
   !if (user%dosoftening) then
   !   if (user%doangle) then 
   !      crater%sinimpang = 0.5_DP * SQRT2 ! Use 45 degree impact angle if user is allowing impact angle to vary
   !   else
   !      crater%sinimpang = 1.0_DP ! Use 90 degree impact angle if user is not allowing impact angle to vary
   !   end if
   !   crater%impvel = rmsvel
   !   do k=2,domain%pnum
   !      crater%imp = prod(1,k)
   !      dN = (prod(2,k - 1) - prod(2,k)) / domain%area / user%interval
   !      ! Do it as the above code trying to figure out the size of a crater with both strength models
   !      crater%strflag = 0 
   !      call crater_generate(user,crater,domain)
   !      if (crater%fcrat > domain%smallest_crater) exit
   !      ! First get the baseline per-crater diffusion rate
   !      kappaNdot = 0.1_DP * crater%fcrat**4
!
!         ! Now add in the extra per-crater diffusion
!         kappaNdot = kappaNdot + (SOFTEN_FACTOR * 0.25_DP * PI) * crater%fcrat**(2._DP + SOFTEN_SLOPE)
!         domain%subpixel_diffusion_const = domain%subpixel_diffusion_const + kappaNdot * dN
!      end do
!   end if
   domain%subpixel_diffusion_const = 0.0_DP
   !write(*,*) "Sub-pixel kappa: ",domain%subpixel_diffusion_const

   domain%initialize = .false.
   return

end subroutine init_domain
