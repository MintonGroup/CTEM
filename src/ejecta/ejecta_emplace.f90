!****f* ejecta/ejecta_emplace
! Name
!   ejecta_emplace -- Calculate ejecta mass during excavation stage.
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_io
!   * module_crater
!   * module_regolith
!   * module_ejecta
!   
!   call ejecta_emplace(user,surf,crater,domain,ejb,ejtble)
!
! DESCRIPTION
!    
!   The estimation of ejecta mass is to treat a crater cavity as parabola shell (Richardson 2009). 
!   At a given landing distance of an ejecta block, one can always find its origin within a transient
!   crater, and as a result, one can know the shape of its parabola shell. 
!
!   It also includes stream tube's calculation. Yet, we use ejecta mass that is obtained from 
!   parabola shell approximation as a constraint of the stream tube's shape. 
! 
!   Besides, we improved ejecta mass distribution by adopting a ray model based on Superformula. 
!   Citation: Gielis, J. "A Generic Geometric Transformation that Unifies a Wide Range of Natural 
!   and Abstract Shapes." Amer. J. Botany 90, 333-338, 2003. 
! 
! ARGUMENTS
!   Input
!   * user    -- User input parameters
!   * surf    -- Surface ggrid
!   * crater  -- Crater dimension container
!   * domain  -- Simulation domain variable container
!   * ejb     -- Ejecta blanket lookup table
!   * ejtble  -- Ejecta blanket lookup table length 
!
!   Output
!   * surf    -- Outputs the new ejecta blanket onto the grid
!   * crater  -- May affects the value of the maximum affected distance
!   * ejbmass -- Outputs total ejecta thickness for mass conservation in crater_emplace
! 
! Notes
!   The cutoff value of ejecta mass for the stream tube volume's calculation is hard coded. 
!
!***

!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_emplace
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Emplaces the ejecta around the crater and softens the terrain using the topographic diffusion model (if the user
!  requests). 
!  
!
!  Input
!    Arguments : user   :     The user-defined variables from the input file
!                surf   :     Surface grid 
!                crater :     Crater dimension container
!                domain :     Simulation domain variables container
!                ejb    :     Ejecta blanket lookup table 
!                ejtble :     Ejecta blanket lookup table length
!
!  Output
!    Arguments : surf   :     Outputs the new ejecta blanket onto the grid
!                crater :     May affects the value of the maximum affected distance
!           
! 
!  Notes       : Crater ray model is based on a mathematic formula, superformula, developed by a belgium scientist, Johan Gielis. It
!                attempts to simulate the shapes of biological creatures such as sea animals (starfish) or bacteria. This finding was 
!                getting attention from Nature and Science magzine. Citation: Gielis, J. "A Generic Geometric Transformation that Unifies
!                a Wide Range of Natural and Abstract Shapes." Amer. J. Botany 90, 333-338, 2003.
!                The cutoff of ejecta thickness is still buggy.  
!
!**********************************************************************************************************************************
subroutine ejecta_emplace(user,surf,crater,domain,ejb,ejtble,ejbmass)
   use module_globals
   use module_util
   use module_io
   use module_crater
   use module_regolith
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_emplace
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   integer(I4B),intent(in) :: ejtble
   type(ejbtype),dimension(ejtble),intent(in)    :: ejb
   real(DP),intent(out) :: ejbmass

   ! Internal variables
   real(DP) :: lrad,lrad0,lradsq,cdepth,distance,erad,craterslope,baseline,inneredge,outeredge
   integer(I4B),parameter :: MAXLOOP = 10 ! Maximum number of times to loop the ejecta angle correction calculation
   integer(I4B) :: xpi,ypi,i,j,k,n,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq,radsq,ebh,ejdissq,continuous
   real(DP),dimension(:,:),allocatable :: cumulative_elchange,big_cumulative_elchange
   integer(I4B),dimension(:,:,:),allocatable :: indarray,big_indarray
   integer(I4B) :: bigi,bigj
   character(len=MESSAGESIZE) :: message  ! message for the progress bar
   !Ray model parameters and variables by Gielis superformula
   integer(I4B) :: nrays
   real(DP),dimension(15) :: rn
   real(DP) :: theta, lradp
   real(DP), parameter :: n1 = 4.0_DP
   real(DP) :: n2, mag
   ! Ray Mass conservation
   real(DP), parameter :: a = 16.8799 !a = 11.8126 ! Fitting parameters for a relation between ray length and radius of crater 
   real(DP), parameter :: b = 0.120621 !0.143   ! based on Jake's crater rays mapping studies! 
   real(DP) :: mvrld                  ! median value of ray length distribution
   real(DP) :: mvrldsc                ! median value of ray length distribution scaled by continuous ejecta extent

   ! Enhanced factor test
   real(DP)     :: xef, yef, thetamax
   integer(I4B) :: xefpi, yefpi, nef, ray_pix
   integer(I4B), dimension(:), allocatable :: sf, tot
   real(DP), dimension(:), allocatable     :: ef
 
   ! Flowery ray variables
   integer(I4B)  :: nfrays
   real(DP)      :: n1f, n2f, magf, lradf, rayf
   ! Ray mixing model parameters:
   real(DP)      :: h_raymix ! Ray mixing depth: l = 0.5755 * (D_sc)^-0.3136 * R_p^1.25, D_sc = 8 * h
                             !                   h = 0.021 * l^-3.188 * R_p^3.985 in unit of kilometers
   real(DP)      :: h_raymixratio ! Ray mixing ratio: h / L = h / (0.1 R_p^0.74 * (l/R_p)^-4.37
                                  ! h/L = 0.2137 * R_p^0.052 * (l/R_p)^1.182
   real(DP), parameter :: k_raymixratio = 0.171726_DP !1.7172_DP 
   real(DP), parameter :: b_lrad1  = -3.188_DP !1.181_DP
   real(DP), parameter :: b_frad1  =  0.79719_DP !0.057_DP
   real(DP), parameter :: SCD = 0.125_DP
   real(DP) :: vsq, ejtheta, melt, vol_sc, dsc

   ! Melt zone's radius
   real(DP) :: rm, dm
   !call regolith_melt_zone(user,crater,crater%imp,crater%impvel,rm,dm)
   !write(*,*) crater%imp, crater%frad, rm
   
   ! Executable code

   cdepth = DDRATIO * crater%fcrat
   crater%vdepth = crater%ejrim + cdepth
   crater%vrim   = crater%ejrim + crater%rheight
   
   if (crater%ejdis <= crater%rad) return

   ! determine area to effect
   inc = max(min(crater%ejdispx + 1,PBCLIM*user%gridsize),1)
   crater%maxinc = max(crater%maxinc,inc)
   radsq = crater%rad**2
   incsq = inc**2
   ejdissq = crater%ejdis**2

   if (inc >= user%gridsize / 2) then
      if (user%testflag) then
          write(*,*) 'Big ejecta: fcrat =',crater%fcrat, ' Ej/S =',(crater%ejdispx*user%pix)/domain%side, ' Ejrim =', crater%ejrim
       else
         write(message,'("Ejb: Dc=",ES9.2," Ej/S=",F0.3)') crater%fcrat,(crater%ejdispx*user%pix)/domain%side
         call io_updatePbar(message)
       end if
   endif
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   allocate(indarray(2,-inc:inc,-inc:inc))

   call random_number(rn)
   cumulative_elchange = 0._DP
   indarray = inc ! initialize this array to point to a corner (this should have 0 elevation change since we're only doing work
                ! within a circle of radius irad

   ! *************************** Continuous Ejecta Formula  *****************************!
   continuous = 2.3_DP * crater%frad**(1.006_DP)

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
   ! There are three parameters regarding to our desired ray shape: n1, n2, m
   ! m:  the repeating part of formula, and it controls the number of rays (arms).
   ! n1: the default value is set as 4.0 in our rays case
   ! n2: n2 and n1 combining together is to control the slenderness of a ray, in general, n1/n2 is always smaller than 1 in our cases.
   !     The smaller the ratio, the skinnier a ray.

   nrays      = nint(8 * rn(15)) + 6
   mvrld      = crater%ejdis
   mvrldsc    = mvrld / continuous
   n2         = 8.0_DP * ( log10(mvrldsc) / log10(2.0_DP) ) + 2.0_DP

   ! *************************** Part II. Flowery Ray Model         ************************************!
   nfrays = 20
   n1f  = 1.0_DP
   n2f  = 0.5_DP
   rayf = 7.5_DP
 
   ! Enhanced factor test: Build a enhanced factor quick lookup table
   ! Use a circular sector with an angle containing a half ray, which is pi/nrays. First, we determine the looping area for 
   ! this sector, then the look-up table will be built while looping over each pixel and then making it to go to the right bin. 
   ! The size of bin is one pixel. 
   xef        = mvrld   
   thetamax   = PI/dble(nrays)                   
   yef        = mvrld * sin(thetamax)
   xefpi      = nint(xef / user%pix)
   yefpi      = nint(yef / user%pix)
   nef        = ceiling( (mvrld - crater%rad) / user%pix)
   allocate(sf(nef))
   allocate(tot(nef))
   allocate(ef(nef))
   tot = 0
   sf  = 0

   do j = yefpi, 0, -1
      do i = xefpi, 0, -1
         xp     = (i + crater%xlpx) * user%pix
         yp     = (j + crater%ylpx) * user%pix
         lradsq = (xp - crater%xl)**2 + (yp - crater%yl)**2
         lrad   = sqrt(lradsq) 
         ! inside or outside a ray?
         theta  = atan2(j * 1._DP,i * 1._DP)
         mag    = ( (abs(cos(nrays * theta / 4.0_DP)))**n2 + (abs(sin(nrays * theta / 4.0_DP)))**n2 )**(-1.0_DP/n1) !&
         
         magf   = rayf * ( (abs(cos(nfrays * theta / 4.0_DP)))**n2f + (abs(sin(nfrays * theta / 4.0_DP)))**n2f )**(-1.0_DP/n1f) 
         lradp  = continuous * mag
         lradf  = continuous * magf 
         lradp  = max(lradp, lradf)
         !if ( lrad > continuous .and. lrad < mvrld .and. theta < thetamax .and. theta>0._DP) then
         if (lrad >= crater%rad .and. lrad < mvrld .and. theta < thetamax .and. theta>0._DP) then
            ray_pix  = ceiling( (lrad - crater%rad) / user%pix)
            tot(ray_pix) = tot(ray_pix) + 1
            if (lrad < lradp) sf(ray_pix) = sf(ray_pix) + 1
         end if
      end do
   end do

   ! Smooth out the enhanced factor lookup table
   do i=1,nef
      ef(i) = max( dble(tot(i)) / dble(sf(i)), 1.0_DP)
   end do

   do i=1,nef
      if (ef(i) /= ef(i) .or. ef(i) > VBIG) then 
         ef(i) = dble(tot(i)) / 1.0 
      end if
   end do
   
   deallocate(sf)
   deallocate(tot)
   outeredge = crater%frad + domain%ejbres * (EJBTABSIZE - 0.5_DP)
   inneredge = crater%frad + 0.5_DP * domain%ejbres 
  
   ejbmass = 0.0_DP
!   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
!   !$OMP SHARED(user,domain,crater,surf,ejb,ejtble,mvrld,mvrldsc,nrays,n2,nef,ef,n2f,n1f,nfrays,rayf) &
!   !$OMP SHARED(inc,incsq,ejdissq,fradsq,radsq,indarray,cumulative_elchange,rn,continuous) 
!   !$OMP REDUCTION(+:massray,massej,massrayef)
   do j = -inc,inc
      do i = -inc,inc

         ! find distance from crater center
         iradsq = i*i + j*j
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! Find distance from crater center to current pixel center in real space
         xp = xpi * user%pix
         yp = ypi * user%pix
         
         !lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2
         !lrad = sqrt(lradsq)

            

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)

         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi

         ! Find angle-corrected landing distance
         distance = sqrt((crater%xl - xp)**2 + (crater%yl - yp)**2)
         baseline = crater%melev + ((i * crater%xslp) + (j * crater%yslp)) * user%pix - surf(xpi,ypi)%dem
         craterslope = atan(baseline / distance)
         lrad = distance
         lrad0 = lrad
         do n = 1,MAXLOOP
            k = max(min(1 + int((lrad - inneredge) / (outeredge - inneredge) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
            call ejecta_interpolate(crater,domain,lrad,ejb,ejtble,ebh,vsq=vsq,theta=ejtheta,erad=erad)
            lrad =  distance * (erad + sqrt(vsq) * sin(2*ejtheta)) / (erad + sqrt(vsq) * (sin(2 * (ejtheta + craterslope))))
            if (abs(lrad - lrad0) < domain%small) exit
            lrad0 = lrad
         end do
         lradsq = lrad**2


         if ((lradsq <= ejdissq) .and. (lradsq >= radsq) .and. (lrad > 0.0_DP)) then
            theta = atan2(j * 1._DP,i * 1._DP) + 2.0_DP * PI
            mag   = ( ( (abs(cos(nrays * theta / 4.0_DP)))**n2 + &
                    (abs(sin(nrays * theta / 4.0_DP)))**n2 )**(-1.0_DP/n1) ) 
            magf  = rayf * ( ( (abs(cos(nfrays * theta / 4.0_DP)))**n2f + &
                      (abs(sin(nfrays * theta / 4.0_DP)))**n2f )**(-1.0_DP/n1f))
            lradp = continuous * mag
            lradf = continuous * magf
            lradp = max(lradp, lradf) 
               

            if (lrad < lradp) then
               call ejecta_interpolate(crater,domain,lrad,ejb,ejtble,ebh,vsq=vsq,theta=ejtheta,melt=melt)
               ray_pix = ceiling((lrad - crater%rad) / user%pix)
               ebh       = ebh * ef(ray_pix)

               if (user%doregotrack .and. ebh>1.0e-8) then
                  call regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,rm)
                  !call regolith_transport(user,surf(xpi,ypi),crater,domain,ejb,ejtble,lrad,ebh,newlayer)
                  !print *, lrad / crater%frad, comp 
                  dsc = ebh + SCD * 1.161_DP * (ebh**0.78) * (sqrt(vsq)**0.44) * (user%gaccel**(-0.22)) * (sin(ejtheta)**(1.0/3.0))
                  if (dsc - ebh > 1.0e-08) then
                     call regolith_mix(surf(xpi,ypi), dsc)
                  end if
                  !vol_sc = PI / 48.0 * dsc**3
                  !h_raymixratio = dsc / ebh
                  !if (i > 0 .and. j > 0 .and. lrad > continuous) write(*,*) lrad/crater%frad, ebh, sqrt(vsq), &
                  !dsc, h_raymixratio  
                  ! Check out: Ray_Mixing_Model.dox
                  !h_raymixratio = ALPHA * k_raymixratio * (lrad/crater%frad)**(b_lrad1) * (crater%frad/1000.0_DP)**(b_frad1)
                  !h_raymix = h_raymixratio * ebh 
                  !if (h_raymix > ebh) then
                  !   call regolith_mix(surf(xpi,ypi), h_raymix)
                  !end if
               end if

            else
               ebh = 0._DP
            end if

            cumulative_elchange(i,j) = cumulative_elchange(i,j) + ebh
            ejbmass = ejbmass + ebh
         end if
         
      end do
   end do
!   !$OMP END PARALLEL DO

   deallocate(ef)

   ! Create box for soften calculation (will be no bigger than the grid itself)
   if (2 * inc + 1 < user%gridsize) then
      call ejecta_soften(user,surf,2 * inc + 1,indarray,cumulative_elchange)
      ! Add the ejecta back to the DEM
      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j)
            ypi = indarray(2,i,j)
            surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cumulative_elchange(i,j) 
            surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cumulative_elchange(i,j), 0.0_DP)
         end do
      end do
   else ! Ejecta wraps around the grid. 
        ! We will therefore send in the whole grid with the total ejecta thickness added to each pixel
      allocate(big_cumulative_elchange(0:user%gridsize+1,0:user%gridsize+1))
      allocate(big_indarray(2,0:user%gridsize+1,0:user%gridsize+1))
      do bigj = 0,user%gridsize + 1
         do bigi = 0,user%gridsize + 1
            xpi = bigi
            ypi = bigj
            call util_periodic(xpi,ypi,user%gridsize)
            big_indarray(1,bigi,bigj) = xpi
            big_indarray(2,bigi,bigj) = ypi
            big_cumulative_elchange(bigi,bigj) = 0.0_DP
         end do
      end do

      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j) 
            ypi = indarray(2,i,j)
            big_cumulative_elchange(xpi,ypi) = big_cumulative_elchange(xpi,ypi) + cumulative_elchange(i,j)
         end do
      end do
      call ejecta_soften(user,surf,user%gridsize + 2,big_indarray,big_cumulative_elchange)
      do i = 1,user%gridsize
         do j = 1,user%gridsize
            surf(i,j)%dem = surf(i,j)%dem + big_cumulative_elchange(i,j) 
            surf(i,j)%ejcov = max(surf(i,j)%ejcov + big_cumulative_elchange(i,j),0.0_DP)
         end do
      end do
      deallocate(big_cumulative_elchange,big_indarray)
   end if

   !if (user%doregotrack) call regolith_rays(user,crater,domain,ejtble,ejb)

   deallocate(cumulative_elchange,indarray)

   return
end subroutine ejecta_emplace

