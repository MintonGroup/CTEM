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
subroutine ejecta_emplace(user,surf,crater,domain,ejb,ejtble,deltaMtot,age,age_resolution)
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
   real(DP),intent(in) :: deltaMtot
   real(DP),intent(in)  :: age
   real(DP),intent(in)  :: age_resolution

   ! Internal variables
   real(DP) :: lrad,lradsq,cdepth
   integer(I4B),parameter :: MAXLOOP = 100 ! Maximum number of times to loop the ejecta angle correction calculation
   integer(I4B) :: xpi,ypi,i,j,k,n,inc,incsq,iradsq,idistorted,jdistorted
   real(DP) :: xp,yp,fradsq,fradpxsq,radsq,ebh,ejdissq,ejbmass,fmasscons,areafrac,xbar,ybar,krad,kdiffmax
   real(DP),dimension(:,:),allocatable :: cumulative_elchange,big_cumulative_elchange,kdiff,big_kdiff,cel,big_cel
   integer(I4B),dimension(:,:,:),allocatable :: indarray,big_indarray
   real(DP),dimension(:,:),allocatable :: ejdistribution,diffdistribution
   integer(I4B) :: bigi,bigj,maxhits,nin,nnot,dradsq
   character(len=MESSAGESIZE) :: message  ! message for the progress bar
   

   ! Ray mixing model variables 
   real(DP)      :: dsc

   ! Melt zone's radius
   real(DP) :: rm, dm, melt, eradc

   ! Ejecta pattern distortion parameters
   real(DP) :: distance,erad,craterslope,landslope,baseline,lrange,frac,ejheight,ebh0,maxdistance
   real(DP) :: maxslp
   real(DP)      :: vsq, ejtheta
   integer(I4B) :: ind,klo

   ! Age
   real(SP) :: age_mean


   ! Executable code

   call regolith_melt_zone(user,crater,crater%imp,crater%impvel,rm,dm)

   cdepth = DDRATIO * crater%fcrat
   crater%vdepth = crater%ejrim + cdepth
   crater%vrim   = crater%ejrim + crater%rheight
   
   if (crater%ejdis <= crater%rad) return

   ! determine area to effect
   inc = max(min(nint(min(PI * user%trad / user%pix, min(crater%ejdis, crater%frad * user%ejecta_truncation / user%pix))) + 1, &
            user%gridsize - nint(crater%continuous / user%pix)),1)
   maxdistance = inc * user%pix

   ! Increase the box a bit to take into account possible ejecta pattern distortion due to topography
   inc = ceiling(inc * 1.5_DP)

   if (user%dosoftening) then
      krad = user%ejecta_truncation * crater%frad
      kdiffmax = user%Kd1 * crater%frad**(user%psi)
      dradsq = int(krad / user%pix) + 3
      inc = max(inc,dradsq)
      dradsq = dradsq**2
   end if

   crater%maxinc = max(crater%maxinc,inc)
   radsq = crater%rad**2
   fradsq = crater%frad**2
   fradpxsq = crater%fradpx**2
   ejdissq = crater%ejdis**2

   incsq = inc**2

   if (inc >= user%gridsize / 2) then
      if (user%testflag) then
          write(*,*) 'Big ejecta: fcrat =',crater%fcrat, ' Ej/S =',(crater%ejdispx*user%pix)/domain%side, ' Ejrim =', crater%ejrim
       else
         write(message,'("Ejb: Dc=",ES9.2," Ej/S=",F0.3)') crater%fcrat,(crater%ejdispx*user%pix)/domain%side
         call io_updatePbar(message)
       end if
   endif
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   allocate(cel(-inc:inc,-inc:inc))
   allocate(kdiff(-inc:inc,-inc:inc))
   allocate(indarray(2,-inc:inc,-inc:inc))
   allocate(ejdistribution(-inc:inc,-inc:inc))
   allocate(diffdistribution(-inc:inc,-inc:inc))

   cumulative_elchange = 0.0_DP
   kdiff = 0.0_DP
   indarray = inc - 1 ! initialize this array to point to a corner (this should have 0 elevation change since we're only doing work
                ! within a circle of radius irad

   ! *************************** Continuous Ejecta Formula  *****************************!
   call ejecta_ray_pattern(user,surf,crater,inc,-inc,inc,-inc,inc,diffdistribution,ejdistribution)

   ejbmass = 0.0_DP
   nin = 0
   nnot = 0

   !!$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !!$OMP SHARED(user,domain,crater,surf,ejb,ejtble) &
   !!$OMP SHARED(inc,incsq) &
   !!$OMP SHARED(cumulative_elchange,kdiff,kdiffmax,indarray,ejdistribution,diffdistribution) 
   do j = -inc,inc
      do i = -inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! Find distance from crater center to current pixel center in real space
         xp = xpi * user%pix
         yp = ypi * user%pix

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)

         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi

         lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2
         lrad = sqrt(lradsq)

         ! Estimate ejecta pattern distortion due to target surface angle and topography
         ! This must be done iteratively because the ejection distance and ejection angle vary
         distance = lrad
         maxslp = -huge(maxslp)
         klo = int((log(lrad) - log(crater%rad)) / domain%ejbres)
         do n = 1,MAXLOOP
            call ejecta_interpolate(crater,domain,distance,ejb,ejtble,ebh,vsq=vsq,theta=ejtheta,erad=erad,melt=melt)
            if ((n > 1).and.((abs(ebh0 - ebh) / ebh0) < domain%small)) exit
            ebh0 = ebh
               
            erad = exp(erad)
            lrange = lrad - erad

            baseline = ((i * crater%xslp) + (j * crater%yslp)) * user%pix
            craterslope = atan(baseline / lrad)
            if (craterslope > maxslp) maxslp = craterslope

            ejheight = erad * sin(craterslope) + crater%melev
           
            landslope = atan((surf(xpi,ypi)%dem - ejheight) / lrange)

            ! Calculate corrected landing velocity for this location
            vsq = (lrange * user%gaccel * cos(craterslope)) / &
                  (sin(2 * ejtheta + 2 * craterslope - landslope) - sin(landslope))
            if (vsq < 0.0_DP) exit

            ! Find out where in the table this new velocity corresponds to
            ind = 1
            call util_search(ejb%vesq,ind,ejtble,vsq,klo)
            klo = min(max(klo,1),ejtble-1)
            ! Interpolate on the table to find the flat plane equivalent landing distance for this velocity
            frac = (vsq - ejb(klo)%vesq) / (ejb(klo+1)%vesq - ejb(klo)%vesq)
            distance = exp(ejb(klo)%lrad) + frac * (exp(ejb(klo+1)%lrad) - exp(ejb(klo)%lrad))
         end do 

         if (vsq < 0.0_DP) cycle
         if (distance /= distance) cycle
         idistorted = int(i * distance / lrad)
         if (abs(idistorted) > inc) cycle
         jdistorted = int(j * distance / lrad)
         if (abs(jdistorted) > inc) cycle
        
         iradsq = idistorted**2 + jdistorted**2
         if ((iradsq > incsq).or.(distance <= crater%rad)) cycle

         ! we need to cut a hole out from the inside of the crater
         xbar = xpi * user%pix - crater%xl 
         ybar = ypi * user%pix - crater%yl

         areafrac =  (1.0_DP - util_area_intersection(crater%rad,xbar,ybar,user%pix)) 

         ebh = areafrac * ejdistribution(idistorted,jdistorted) * ebh
         cumulative_elchange(i,j) = areafrac * cumulative_elchange(i,j) + ebh

         if (user%dosoftening) then
            ! Do extra diffusive degradation over ejecta region
            areafrac =  (1.0_DP - util_area_intersection(crater%frad,xbar,ybar,user%pix)) 
            kdiff(i,j) = areafrac * diffdistribution(idistorted,jdistorted) * kdiffmax
         end if


         if (user%doregotrack .and. ebh>1.0e-8_DP) then
            call regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,rm,vsq,age,age_resolution)
         end if

            
      end do
   end do
   !!$OMP END PARALLEL DO
   ejbmass = sum(cumulative_elchange)

   ! Create buffer to prevent infinite hole bug
   kdiff(-inc,:) = 0.0_DP
   kdiff(inc,:) = 0.0_DP
   kdiff(:,-inc) = 0.0_DP
   kdiff(:,inc) = 0.0_DP

   cumulative_elchange(-inc,:) = 0.0_DP
   cumulative_elchange(inc,:) = 0.0_DP
   cumulative_elchange(:,-inc) = 0.0_DP
   cumulative_elchange(:,inc) = 0.0_DP

   ! Do mass conservation by adjusting ejecta thickness
   fmasscons = (-deltaMtot)/ ejbmass
   cumulative_elchange = cumulative_elchange * fmasscons
   maxhits = 1
   ! Create box for soften calculation (will be no bigger than the grid itself)
   if (2 * inc + 1 < user%gridsize) then

      ! extra soften calculation
      if (user%dosoftening) then
         cel = 0.0_DP
         call util_diffusion_solver(user,surf,2 * inc + 1,indarray,kdiff,cel,maxhits)
         do j = -inc,inc
            do i = -inc,inc
               xpi = indarray(1,i,j)
               ypi = indarray(2,i,j)
               surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cel(i,j)
               surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cel(i,j),0.0_DP)
            end do
         end do
      end if

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
      allocate(big_cel(0:user%gridsize+1,0:user%gridsize+1))
      allocate(big_indarray(2,0:user%gridsize+1,0:user%gridsize+1))
      allocate(big_kdiff(0:user%gridsize+1,0:user%gridsize+1))

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

      big_kdiff = 0.0_DP

      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j) 
            ypi = indarray(2,i,j)
            big_cumulative_elchange(xpi,ypi) = big_cumulative_elchange(xpi,ypi) + cumulative_elchange(i,j)
            big_kdiff(xpi,ypi) = big_kdiff(xpi,ypi) + kdiff(i,j)
         end do
      end do

      if (user%dosoftening) then
         big_cel = 0.0_DP
         call util_diffusion_solver(user,surf,user%gridsize + 2,big_indarray,big_kdiff,big_cel,maxhits)

         do j = 1,user%gridsize
            do i = 1,user%gridsize
               surf(i,j)%dem = surf(i,j)%dem + big_cel(i,j)
               surf(i,j)%ejcov = max(surf(i,j)%ejcov + big_cel(i,j),0.0_DP)
            end do
         end do
      end if

      call ejecta_soften(user,surf,user%gridsize + 2,big_indarray,big_cumulative_elchange)

      do i = 1,user%gridsize
         do j = 1,user%gridsize
            surf(i,j)%dem = surf(i,j)%dem + big_cumulative_elchange(i,j) 
            surf(i,j)%ejcov = max(surf(i,j)%ejcov + big_cumulative_elchange(i,j),0.0_DP)
         end do
      end do

      deallocate(big_cumulative_elchange,big_indarray,big_kdiff,big_cel)
   end if

   !if (user%doregotrack) call regolith_rays(user,crater,domain,ejtble,ejb)

   deallocate(cumulative_elchange,indarray,diffdistribution,ejdistribution,kdiff,cel)

   return
end subroutine ejecta_emplace

