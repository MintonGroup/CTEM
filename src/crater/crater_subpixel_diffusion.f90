!**********************************************************************************************************************************
!
!  Unit Name   : crater_subcrater_diffusion
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
subroutine crater_subpixel_diffusion(user,surf,nflux,domain,finterval,kdiffin)
   use module_globals
   use module_util
   use module_ejecta
   use module_crater, EXCEPT_THIS_ONE => crater_subpixel_diffusion
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   real(DP),dimension(:,:),intent(in) :: nflux 
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: finterval
   real(DP),dimension(:,:),intent(inout) :: kdiffin

   ! Internal variables
   real(DP),dimension(0:user%gridsize + 1,0:user%gridsize + 1) :: cumulative_elchange,kdiff
   integer(I4B),dimension(2,0:user%gridsize + 1,0:user%gridsize + 1) :: indarray
   integer(I4B) :: i,j,l,xpi,ypi,k,ktot,kk,ioerr,inc,incsq,iradsq,xi,xf,yi,yf,crat,imin,imax,jmin,jmax
   integer(I8B) :: m,N
   integer(I4B) :: maxhits = 1
   real(DP) :: dburial,lambda,dKdN,diam,radius,Area,avgejc
   real(DP) :: dN,diam_regolith,diam_bedrock,RR
   real(DP),dimension(3)    :: rn  
   real(DP) :: superlen,rayfrac,cutout,fe,fd
   type(cratertype) :: crater
   real(DP),dimension(:,:),allocatable :: diffdistribution,ejdistribution
   integer(I4B),dimension(:,:),allocatable :: ejisray
   real(DP) :: xbar,ybar,dD,xp,yp,areafrac,krad,supersize,lrad
   integer(I8B),dimension(user%gridsize,user%gridsize) :: Ngrid
   

   ! Create box for soften calculation (will be no bigger than the grid itself)
   do j = 0,user%gridsize + 1
      do i = 0,user%gridsize + 1
         xpi = i
         ypi = j
         call util_periodic(xpi,ypi,user%gridsize)
         indarray(1,i,j) = xpi
         indarray(2,i,j) = ypi
         kdiff(i,j) = kdiffin(xpi,ypi)
      end do
   end do

   avgejc = sum(surf%ejcov) / domain%parea

   fe = FEPROX
   fd = user%ejecta_truncation
   if (user%dosoftening) fe = user%fe

   ! Generate both the subpixel and superdomain diffusive degradation
   do k = 1,domain%pnum - 1
      dN = nflux(3,k) * user%interval * finterval

      diam_bedrock = nflux(1,k)
      diam_regolith = nflux(2,k)

      dburial = 0.5_DP * EXFAC * max(diam_bedrock,diam_regolith)
      if (dburial > avgejc) then
         diam = diam_bedrock
      else
         diam = diam_regolith
      end if
      radius = 0.5_DP * diam

      if ((fd * diam < user%pix) .or. (dN * PI * (fe * radius)**2 > 0.1_DP)) then 
      !Do the average degradation per pixel for the subpixel component
         !Empirically-derived "intrinsic" degradation function from proximal ejecta redistribution
         dKdN = KD1PROX * PI * FEPROX**2 * (radius)**(2.0_DP + PSIPROX) / domain%parea
         if (user%dosoftening) then 
         ! User-defined degradation function
            dKdN = dKdN + PI * user%fe**2 * radius**2 * crater_degradation_function(user,radius) / domain%parea
         end if

         lambda = dN * domain%parea

         ! Don't parallelize the random
         do j = 1,user%gridsize
            do i = 1,user%gridsize
               Ngrid(i,j) = util_poisson(lambda)
            end do
         end do

         !$OMP PARALLEL DO DEFAULT(PRIVATE) &
         !$OMP SHARED(user) &
         !$OMP SHARED(Ngrid,dKdN,kdiff,dN,lambda) 
         do j = 1,user%gridsize
            do i = 1,user%gridsize
               if ((Ngrid(i,j)) == 0) cycle
               if ((Ngrid(i,j) > 0).and.(Ngrid(i,j) == Ngrid(i,j))) then
                  kdiff(i,j) = kdiff(i,j) + dKdN * Ngrid(i,j)
               else ! In case of integer overflow
                  kdiff(i,j) = kdiff(i,j) + dKdN * lambda
               end if
            end do
         end do
         !$OMP END PARALLEL DO

      else if (user%dosoftening) then
      ! Do the degradation as individual circles

         superlen = fd * diam + domain%side
         cutout = 0.0_DP
         crater%continuous = RCONT * radius**(EXPCONT) 
         if (diam > domain%smallest_crater) then
            if (.not.user%superdomain) exit
            ! Superdomain craters 
            cutout = domain%side + crater%continuous
         end if
         Area = superlen**2
         if (diam < domain%biggest_crater) Area = Area - cutout**2

         lambda = dN * Area
         dD = nflux(2,k+1) - nflux(2,k)

         N = util_poisson(lambda)
         do m = 1, N
            call random_number(rn)

            if (diam < domain%smallest_crater) then
               ! Subpixel craters with large degradation region
               crater%xl = real(superlen * (rn(1) - 0.5_DP) + 0.5_DP * domain%side, kind=SP) 
               crater%yl = real(superlen * (rn(2) - 0.5_DP) + 0.5_DP * domain%side, kind=SP) 
            else if (.not.user%superdomain) then
               cycle
            else
               ! Superdomain craters 
               crater%xl = real((superlen - cutout) * (rn(1) - 0.5_DP), kind=SP) 
               if (crater%xl > 0.0_DP) crater%xl = crater%xl + cutout
               crater%yl = real((superlen - cutout) * (rn(2) - 0.5_DP), kind=SP) 
               if (crater%yl > 0.0_DP) crater%yl = crater%yl + cutout
            end if

            crater%xlpx = nint(crater%xl / user%pix)
            crater%ylpx = nint(crater%yl / user%pix)

            crater%frad = 0.5_DP * (diam + dD * rn(3))
            crater%continuous = RCONT * crater%frad**(EXPCONT) 
            krad = fd * crater%frad
          
            dKdN = crater_degradation_function(user,crater%frad)
            inc  = int(krad / user%pix) + 2
            incsq = inc**2

            xpi = max(crater%xlpx - inc,1)
            imin = xpi - crater%xlpx
            xpi = min(crater%xlpx + inc,user%gridsize)
            imax = xpi - crater%xlpx
            ypi = max(crater%ylpx - inc,1)
            jmin = ypi - crater%ylpx
            ypi = min(crater%ylpx + inc,user%gridsize)
            jmax = ypi - crater%ylpx
 
         
            allocate(diffdistribution(imin:imax,jmin:jmax))
            allocate(ejdistribution(imin:imax,jmin:jmax))
            call ejecta_ray_pattern(user,surf,crater,inc,imin,imax,jmin,jmax,diffdistribution,ejdistribution)
            ! Loop over affected matrix area
            !!$OMP PARALLEL DO DEFAULT(SHARED) IF(inc > INCPAR) &
            !!$OMP FIRSTPRIVATE(jmin,jmax,imin,imax) &
            !!$OMP PRIVATE(i, xp, yp, xpi, ypi, lrad)
            do j = jmin,jmax
               do i = imin,imax
                  xpi = crater%xlpx + i
                  ypi = crater%ylpx + j
                  kdiff(xpi,ypi) = kdiff(xpi,ypi) + dKdN * diffdistribution(i,j)
                  !TEMP
                  xp = xpi * user%pix
                  yp = ypi * user%pix
                  lrad = sqrt((xp - crater%xl)**2 + (yp - crater%yl)**2)
                  surf(xpi,ypi)%ejcov = surf(xpi,ypi)%ejcov + ejdistribution(i,j) &
                  * 0.14_DP * crater%frad**(0.74_DP) * (lrad / crater%frad)**(-3)
               end do
            end do
            !!$OMP END PARALLEL DO
            deallocate(diffdistribution,ejdistribution)
         end do
      end if

   end do

   kdiff(0,0) = kdiff(user%gridsize,user%gridsize)
   kdiff(user%gridsize + 1,user%gridsize + 1) = kdiff(1,1)
   kdiff(0,:) = kdiff(user%gridsize,:)
   kdiff(:,0) = kdiff(:,user%gridsize)
   kdiff(user%gridsize + 1,:) = kdiff(1,:)
   kdiff(:,user%gridsize + 1) = kdiff(:,1)  

   ! write(*,*)
   ! write(*,*) 'avgkdiff = ',sum(kdiff) / (user%gridsize + 2)**2 / finterval
   ! write(*,*)
   ! open(unit=55,file='avgkdiff.dat',status='unknown',position='append')
   ! write(55,*) sum(kdiff) / (user%gridsize + 2)**2 / finterval
   ! close(55)


   call util_diffusion_solver(user,surf,user%gridsize + 2,indarray,kdiff,cumulative_elchange,maxhits)
   do j = 1,user%gridsize
      do i = 1,user%gridsize
         surf(i,j)%dem = surf(i,j)%dem + cumulative_elchange(i,j)
         surf(i,j)%ejcov = max(surf(i,j)%ejcov + cumulative_elchange(i,j),0.0_DP)
      end do
   end do
   kdiffin = 0.0_DP

return
end subroutine crater_subpixel_diffusion

