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
subroutine ejecta_emplace(user,surf,crater,domain,ejb,ejtble)
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

   ! Internal variables
   real(DP) :: lrad,lradsq,cdepth
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq,ebh,ejdissq,continuous
   real(DP),dimension(:,:),allocatable :: cumulative_elchange,big_cumulative_elchange
   integer(I4B),dimension(:,:,:),allocatable :: indarray,big_indarray
   integer(I4B) :: bigi,bigj
   character(len=MESSAGESIZE) :: message  ! message for the progress bar
   !Ray model parameters and variables by Gielis superformula
   integer(I4B) :: nrays
   real(DP),dimension(15) :: rn
   real(DP) :: theta, lradp
   real(DP), parameter :: n1 = 4.0_DP
   real(DP) :: n2, mag, x, y
   ! Ray Mass conservation
   real(DP), parameter :: a = 16.8799 !a = 11.8126 ! Fitting parameters for a relation between ray length and radius of crater 
   real(DP), parameter :: b = 0.120621 !0.143   ! based on Jake's crater rays mapping studies! 
   real(DP) :: mvrld                  ! median value of ray length distribution
   real(DP) :: mvrldsc                ! median value of ray length distribution scaled by continunous ejecta extent
   !real(DP) :: thinnest, lrad_thinnest, v_thinnest, theta_thinnest, rad_sec, vsq

   ! Streamtube
   real(DP) :: comp 
   
   ! Executable code

   cdepth = DDRATIO * crater%fcrat
   crater%vdepth = crater%ejrim + cdepth
   crater%vrim   = crater%ejrim + crater%rheight
   
   if (crater%ejdis <= crater%frad) return

   ! determine area to effect
   inc = max(min(crater%ejdispx + 1,PBCLIM*user%gridsize),1)
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
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

   !Mixing
   !surf(crater%xlpx,crater%ylpx)%nmix = 0
   !TESTING
   !continuous = 23_DP * crater%frad**(1.006_DP)
   continuous = 2.3_DP * crater%frad**(1.006_DP)
   ! From fitting Jake's mapping rays data, it is a linear function about the relationship between the median value of ray length 
   ! distribution and the radius of rayed craters (in unit of kilometers)
   ! Also, we need to scale it with the continuous ejecta's extent for ray model
!#   mvrld      = a * (crater%frad/1000.0)**(b)
!#   mvrldsc    = mvrld / (continuous/crater%frad)
   ! It appears that no strong correlation between the number of rays and size of craters.
   ! The average number of rays is about 10. The minimum and maximum number is 6 and 14 respectively.
   nrays      = nint(8 * rn(15)) + 6
   
   ! Determine parameters of Ray model based on Gielis's superformula
   ! There are three parameters regarding to our desired ray shape: n1, n2, m
   ! m:  the repeating part of formula, and it controls the number of rays (arms).
   ! n1: the default value is set as 4.0 in our rays case
   ! n2: n2 and n1 combining together is to control the slenderness of a ray, in general, n1/n2 is always smaller than 1 in our cases.
   !     The smaller the ratio, the skinnier a ray.
!#   n2 = 8.0_DP * ( log10(mvrldsc) / log10(2.0_DP) ) + 2.0_DP

   ! Testing the maximum ray length regardless of any scales of craters (~58 radii)
   mvrld   = crater%ejdis
   mvrldsc = mvrld / continuous
   n2      = 8.0_DP * ( log10(mvrldsc) / log10(2.0_DP) ) + 2.0_DP
!#   write(*,*) mvrld, mvrldsc, n2, nrays
   !write(*,*) continuous / crater%frad, mvrld
   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !$OMP SHARED(user,domain,crater,surf,ejb,ejtble,mvrld,mvrldsc,nrays,n2) &
   !$OMP SHARED(inc,incsq,ejdissq,fradsq,indarray,cumulative_elchange,rn,continuous)
   do j = -inc,inc
      do i = -inc,inc

         ! find distance from crater center
         iradsq = i*i + j*j
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! Find distance from crater center to current pixel center in real space
         xp = xpi * user%pix
         yp = ypi * user%pix
         
         lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2
         lrad = sqrt(lradsq)

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)

         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi
         
         if ((lradsq <= ejdissq) .and. (lradsq >= fradsq)) then
            ! Model the discontinuous ejecta blanket as a "splat"
            ! Create a splat varying with zimuthal angle with the same length of a splat
            theta = atan2(j * 1._DP,i * 1._DP) + 2.0_DP * PI ! Azimuthal angle
            mag   = ( (abs(cos(nrays * theta / 4.0_DP)))**n2 + (abs(sin(nrays * theta / 4.0_DP)))**n2 )**(1.0_DP/n1)
            x     = cos(theta) / mag
            y     = sin(theta) / mag
            lradp  = continuous * sqrt(x**2 + y**2)
            !lradp = mvrld * crater%frad ! Testing Jake's empirical formula with an assumption of homogeneous ejecta extent 
            if (lrad < lradp) then 
!                We are now in a ray!
               call ejecta_interpolate(crater,domain,lrad,ejb,ejtble,ebh)
!#               write(*,*) i * user%pix / crater%frad, j * user%pix / crater%frad
            else
               ebh = 0._DP
            end if
 
            if (user%doregotrack .and. ebh>1.0e-8) then
               call regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,comp)
               call regolith_transport(user,surf(xpi,ypi),crater,domain,ejb,ejtble,lrad,ebh,comp)
            end if

            cumulative_elchange(i,j) = cumulative_elchange(i,j) + ebh
         end if
         
      end do
   end do
   !$OMP END PARALLEL DO

   if (user%dosoftening) then 
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
   else
      ! Add the ejecta back to the DEM
      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j)
            ypi = indarray(2,i,j)
            surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cumulative_elchange(i,j) 
            surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cumulative_elchange(i,j), 0.0_DP)
         end do
      end do
   end if

   !if (user%doregotrack) call regolith_rays(user,crater,domain,ejtble,ejb)

   deallocate(cumulative_elchange,indarray)

   return
end subroutine ejecta_emplace

