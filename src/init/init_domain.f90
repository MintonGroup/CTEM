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
   type(usertype),intent(inout)                    :: user
   type(cratertype),intent(inout)                  :: crater
   type(domaintype),intent(inout)                  :: domain
   real(DP),dimension(:,:),intent(inout)           :: prod,vdist
   real(DP),dimension(:,:),intent(out)             :: crtscl
   real(DP),dimension(:,:),intent(out),allocatable :: pdist
   real(DP),dimension(:,:),intent(out)             :: nflux

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
   domain%smallest_crater =  user%pix
   domain%smallest_ejecta = SMALLESTEJECTA * user%pix
   domain%smallest_counted_crater = SMALLESTCOUNTABLE * user%pix 
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

   ! Preliminary seismic property calculations
   ! user%seisk = THIRD * user%tvel * user%tfrac
   ! user%cohaccel = user%regcoh / user%trho_r

   ! Now we build an idealized production population
   domain%initialize = .true.

   rmsvel = 0._DP
   numvel = 0._DP
   vlo = -1.0_DP
   do k=1,domain%vnum
      if ((vlo < 0.0_DP).and.(vdist(2,k) > 0._DP)) then
         vlo = vdist(1,k)
         domain%vlo = k
      end if
      if (vdist(3,k) < 1.0_DP) domain%vhi = k 
      rmsvel = rmsvel + vdist(2,k) * vdist(1,k)**2
      numvel = numvel + vdist(2,k)
   end do
   domain%vhi = min(domain%vhi + 1,domain%vnum)
   rmsvel = sqrt(rmsvel/numvel) 
   domain%rmsvel = rmsvel

   domain%smallest_impactor_index = 1
   ! Find the smallest impactor that produces a crater at least 1 pixel wide at the maximum possible impact velocity and angle
   do k = 1,domain%pnum
      crater%imp  = prod(1,k)
      crater%impvel = rmsvel
      if (user%doangle) then 
         crater%sinimpang = 0.5_DP * SQRT2 ! Use 45 degree impact angle if user is allowing impact angle to vary
      else
         crater%sinimpang = 1.0_DP ! Use 90 degree impact angle if user is not allowing impact angle to vary
      end if

      ! Calculate crater size based on both strength models
      crater%strflag = 0
      call crater_generate(user,crater,domain)
      regolithfcrat = crater%fcrat

      crater%strflag = 1
      call crater_generate(user,crater,domain)
      bedrockfcrat = crater%fcrat

      ! Build some useful tables
      nflux(1,k) = bedrockfcrat
      nflux(2,k) = regolithfcrat
      crtscl(1,k) = crater%imp
      crtscl(2,k) = bedrockfcrat

      !Find the smallest impactor to consider
      if (domain%smallest_impactor_index == 1) then
         crater%sinimpang = 1.0_DP ! Use 90 degree impact angle if user is not allowing impact angle to vary
         crater%impvel = vdist(1,domain%vnum) 
         crater%strflag = 0
         call crater_generate(user,crater,domain)
         regolithfcrat = crater%fcrat
         crater%strflag = 1
         call crater_generate(user,crater,domain)
         bedrockfcrat = crater%fcrat
         domain%subcrater_limit = min(regolithfcrat,bedrockfcrat)
         fcrat = max(regolithfcrat,bedrockfcrat) 
         if (fcrat >= user%pix) domain%smallest_impactor_index = k
      end if
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
      
      nflux(3,k) = diffnum / domain%area / user%interval

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

   domain%small = user%pix * SMALLFAC

   domain%initialize = .false.
   return

end subroutine init_domain
