!**********************************************************************************************************************************
!
!  Unit Name   : crater_generate
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Generates random crater
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
subroutine crater_generate(user,crater,domain,prod,vdist,surf)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_generate
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(in),optional :: prod,vdist
   type(surftype),dimension(:,:),intent(in),optional :: surf

   ! Internal variables
   real(DP),dimension(5)    :: rn  ! Random number
   real(DP)                 :: nmark,frac,limp
   real(DP)                 :: lnmark,lprod1,lprod1p,lprod2,lprod2p
   real(DP)                 :: dburial,trfin
   integer(I4B)             :: k,khi,klo,Nk

   ! Get all six random numbers we need in one call
   if (.not.domain%initialize) call random_number(rn)

   ! Find crater center position
   if (domain%initialize) then
      crater%xl = real(domain%side * 0.5_DP, kind=SP)
      crater%yl = real(domain%side * 0.5_DP, kind=SP)
   else if (user%testflag) then
      crater%xl = real(domain%side * 0.5_DP + user%testxoffset, kind=SP)
      crater%yl = real(domain%side * 0.5_DP + user%testyoffset, kind=SP)
   else
      crater%xl = real(domain%side * rn(1), kind=SP)
      crater%yl = real(domain%side * rn(2), kind=SP)
   end if

   crater%xlpx = nint(crater%xl / user%pix)
   crater%ylpx = nint(crater%yl / user%pix) 

   ! Make sure it's on the domain
   call util_periodic(crater%xlpx,crater%ylpx,user%gridsize)

   ! Get impactor size from the distribution
   if (.not.domain%initialize) then
      if (user%testflag) then
         crater%imp = user%testimp
      else 
         if (domain%pnum == 1) then
            crater%imp = prod(1,1)
         else! Draw a random impactor from the production SFD
      
            ! generate random impactor 
            !# nmark = prod(2,domain%smallest_impactor_index) * rn(3)
            nmark = prod(2,domain%smallest_ejecta_index) * rn(3)
            ! Make a guess as to where in the SFD the impactor might be. 
            ! This could speed up the searching if the SFD has a lot of elements in it.
            !# Nk = 1 + domain%pnum - domain%smallest_impactor_index 
            !# klo = domain%smallest_impactor_index
            !# khi = domain%pnum
            !# k = klo + int(Nk * log(rn(3)) / prod(4,khi) / prod(4,klo))
            Nk    = 1 + domain%pnum - domain%smallest_ejecta_index
            klo   = domain%smallest_ejecta_index
            khi   = domain%pnum
            k     = klo + int(Nk * log(rn(3)) / prod(4,khi) / prod(4,klo))
            
            ! Now search the table to find where the impactor actually is
            call util_search(prod,2,domain%pnum-1,nmark,k)
            if (k >= domain%pnum) then
               crater%imp = prod(1,domain%pnum)
            else
               !#if (k <= domain%smallest_impactor_index) k = domain%smallest_impactor_index
               if (k <= domain%smallest_ejecta_index) k = domain%smallest_ejecta_index
               lnmark = log(nmark)
               lprod1 = prod(3,k) !log(prod(1,k))
               lprod1p = prod(3,k + 1) !log(prod(1,k + 1))
               lprod2 = prod(4,k) !log(prod(2,k))
               lprod2p = prod(4,k + 1) !log(prod(2,k + 1))
               frac = (lprod2 - lnmark) / (lprod2 - lprod2p)
               limp = lprod1 + (frac*(lprod1p - lprod1))
               crater%imp = exp(limp)
            end if
         end if
      end if
      crater%imp = crater%imp * (1._DP + 1.0e-3_DP*rn(3)) ! Some user-input SFDs can result in many craters having identical 
                                                          ! diameters. This random number prevents more than one crater from having 
                                                          ! exactly the same diameter, as diameter is used as identification.
   end if
                                                       

   crater%impmass = 4*THIRD*PI*user%prho*(crater%imp*0.5_DP)**3

   ! Determine which strength model to use


   !  find impact angle
   if (.not.domain%initialize) then
      if (user%testflag) then
         crater%sinimpang = sin(user%testang * DEG2RAD)
      else
         if (user%doangle) then
            crater%sinimpang = sqrt(rn(4))
         else
            crater%sinimpang = 1._DP ! Vertical impact only
         end if
      end if
   end if

   if (.not.domain%initialize) then
      if (user%testflag) then
         crater%impvel = user%testvel
      else
         if (domain%vnum == 1) then
            crater%impvel = vdist(1,1)
         else 
            !  Draw impact velocity from the velocity distribution
            nmark = vdist(3,domain%vlo) + (rn(5) * (vdist(3,domain%vhi) - vdist(3,domain%vlo)))
            ! Make a guess as to where in the velocity distribution the impactor might be. 
            ! This could speed up the searching if the velocity distribution has a lot of elements in it.
            klo = domain%vlo
            khi = domain%vhi
            Nk = 1 + khi - klo
            k = domain%vlo + int(Nk * rn(5) * (vdist(3,khi) - vdist(3,klo)))
            call util_search(vdist,3,domain%vnum-1,nmark,k)
            if (k == 0) then
               crater%impvel = vdist(1,1)
            else if (k == domain%vnum) then
               crater%impvel = vdist(1,domain%vnum)
            else
               frac = (nmark-vdist(3,k))/(vdist(3,k+1)-vdist(3,k))
               crater%impvel = vdist(1,k) + (frac*(vdist(1,k+1)-vdist(1,k)))
            end if
         end if
      end if
   end if

   !  scale to crater size
   if (.not.domain%initialize) crater%strflag = 0 ! Begin with regolith strength
   call crater_scale(user,crater%imp,crater%rad,crater%grad,crater%strflag,crater%sinimpang,crater%impvel)
   if (.not.domain%initialize) then
      dburial = EXFAC * crater%rad
      if (dburial > surf(crater%xlpx,crater%ylpx)%ejcov) then
         crater%strflag = 1 ! Use bedrock strength
         call crater_scale(user,crater%imp,crater%rad,crater%grad,crater%strflag,crater%sinimpang,crater%impvel)
      end if
   end if


   trfin = 2 * TRSIM * crater%rad
   if (trfin <= crater%cxtran) then
      crater%fcrat = trfin   ! Simple Crater
   else  ! Complex crater
      crater%fcrat = trfin * ((trfin/crater%cxtran)**crater%cxexp) ! Complex Crater
   end if
   if (crater%imp >= user%basinimp) then
      ! This section is temporary until a better basin scaling law can be implemented. Potter et al. (2012) GRL v39. p 18203
      if ((crater%xl < (0.5_DP * domain%side) ).or.(user%testflag).or.(domain%initialize)) then ! pick TP1 for left-hand hemisphere
         crater%fcrat = 0.1354e3_DP * (0.5_DP * trfin * 1e-3_DP)**(1.389_DP)  ! TP2
         crater%fcrat = 2 * 1.56_DP * crater%fcrat
      else
         crater%fcrat = 0.0718e3_DP * (0.5_DP * trfin * 1e-3_DP)**(1.613_DP)  ! TP1
         crater%fcrat = 2 * 1.2_DP * crater%fcrat
      end if
   end if

   crater%frad = 0.5_DP  * crater%fcrat
   ! Get pixel space values
   crater%fcratpx = nint(crater%fcrat / user%pix)
   return
end subroutine crater_generate

