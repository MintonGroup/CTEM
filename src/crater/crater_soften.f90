!**********************************************************************************************************************************
!
!  Unit Name   : crater_soften
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Softens the terrain of the crater based on empirical studies of equilibrium (see Minton & Fassett 2016)
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
subroutine crater_soften(user,surf,crater,domain)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_soften
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain

   ! Internal variables
   real(DP),dimension(:,:),allocatable :: cumulative_elchange,kdiff
   integer(I4B),dimension(:,:,:),allocatable :: indarray
   real(DP),dimension(0:user%gridsize+1,0:user%gridsize+1) :: big_cel,big_kdiff
   integer(I4B),dimension(2,0:user%gridsize+1,0:user%gridsize+1) :: big_indarray
   integer(I4B) :: maxhits 

   integer(I4B) :: inc,incsq,N,xpi,ypi,iradsq,i,j,ss,ioerr,bigi,bigj
   real(DP) :: xp,yp,fradsq,xbar,ybar,areafrac,kdiffmax,sf

   real(DP) :: mfe,bfe

   ! TEST VARIABLE FE MODEL
   mfe = (user%femax - user%femin) / (log10(user%rmaxfe) - log10(user%rminfe))
   bfe = 0.5_DP * ((user%femax + user%femin) - mfe * (log10(user%rmaxfe) + log10(user%rminfe)))

   crater%fe = mfe * log10(crater%frad) + bfe
   crater%fe = min(user%femax,crater%fe) 
   !^^^^^^^^^^^^^^^^^^^^^^^^^^

   kdiffmax = user%Kd1 * crater%frad**user%psi
   inc = int(min(crater%fe * crater%frad / user%pix,SQRT2 * user%gridsize)) + 3 
   maxhits = (1 + (inc / (user%gridsize / 2)))**2
   crater%maxinc = max(crater%maxinc,inc)
   incsq = (inc - 2)**2

   allocate(kdiff(-inc:inc,-inc:inc))
   allocate(indarray(2,-inc:inc,-inc:inc))
   allocate(cumulative_elchange(-inc:inc,-inc:inc))
   cumulative_elchange = 0.0_DP
   indarray = inc
   kdiff = 0.0_DP

   ! Loop over affected matrix area
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find elevation and grid point
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
            
         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi

         ! find distance from crater center
         iradsq = i*i + j*j

         ! We set the diffusion constant to be proportional to the fraction of pixel area covered by the interior of the crater
         !if (iradsq <= incsq) then
            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
            
            xbar = xp - crater%xl 
            ybar = yp - crater%yl

            areafrac = util_area_intersection(crater%fe * crater%frad,xbar,ybar,user%pix)
             
            kdiff(i,j) = kdiffmax * areafrac  ! Apply user-defined diffusive degradation to degradation region
         !end if
      end do
   end do !end area loopover 
   !write(*,*) 'crater_soften: ',crater%frad,maxval(kdiff)

   if (2 * inc + 1 < user%gridsize) then
      write(*,*) 'crater inc = ',inc
      write(*,*) 'crater kdiff '
      write(*,*) maxval(kdiff)
      call util_diffusion_solver(user,surf,2 * inc + 1,indarray,kdiff,cumulative_elchange,maxhits)
      write(16,*) crater%frad
      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j)
            ypi = indarray(2,i,j)
            write(16,*) i,j,xpi,ypi,kdiff(i,j),surf(xpi,ypi)%dem,cumulative_elchange(i,j)
            !surf(xpi,ypi)%dem = surf(xpi,ypi)%dem + cumulative_elchange(i,j)
            !surf(xpi,ypi)%ejcov = max(surf(xpi,ypi)%ejcov + cumulative_elchange(i,j),0.0_DP)
         end do
      end do
   else
      maxhits = 1
      big_kdiff = 0.0_DP
      big_cel = 0.0_DP
      do bigj = 0,user%gridsize + 1
         do bigi = 0,user%gridsize + 1
            xpi = bigi
            ypi = bigj
            call util_periodic(xpi,ypi,user%gridsize)
            big_indarray(1,bigi,bigj) = xpi
            big_indarray(2,bigi,bigj) = ypi
         end do
      end do

      do j = -inc,inc
         do i = -inc,inc
            xpi = indarray(1,i,j) 
            ypi = indarray(2,i,j)
            big_kdiff(xpi,ypi) = big_kdiff(xpi,ypi) + kdiff(i,j)
         end do
      end do

      big_kdiff(0,0) = big_kdiff(user%gridsize,user%gridsize)
      big_kdiff(user%gridsize + 1,user%gridsize + 1) = big_kdiff(1,1)
      big_kdiff(0,:) = big_kdiff(user%gridsize,:)
      big_kdiff(:,0) = big_kdiff(:,user%gridsize)
      big_kdiff(user%gridsize + 1,:) = big_kdiff(1,:)
      big_kdiff(:,user%gridsize + 1) = big_kdiff(:,1)

      call util_diffusion_solver(user,surf,user%gridsize + 2,big_indarray,big_kdiff,big_cel,maxhits)
      do j = 1,user%gridsize
         do i = 1,user%gridsize
            surf(i,j)%dem = surf(i,j)%dem + big_cel(i,j)
            surf(i,j)%ejcov = max(surf(i,j)%ejcov + big_cel(i,j),0.0_DP)
         end do
      end do
   end if

   deallocate(kdiff,indarray,cumulative_elchange)
return
end subroutine crater_soften

