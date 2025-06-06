!****f* regolith/regolith_streamtube
! Name
!   regolith_streamtube -- Calculate stream tube's volume during excavation stage.
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_regolith
!   
!   call regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,rm)
!
! DESCRIPTION
!    
!   The stream tube is based on Maxwell Z-model, and we were able to derive an analytical function 
!   for a three dimensional stream tube. As a result, we can analyze segments' properties that is 
!   included in a stream tube during excavation stage.
!
! ARGUMENTS
!   Input
!   * user     -- The user-defined variables from the input file
!   * surf     -- Surface grid
!   * crater   -- Crater dimension container
!   * domain   -- Simulation domain variable container
!   * ejb      -- Ejecta blanket lookup table
!   * ejtble   -- Ejecta blanket lookup table length
!   * xp,yp    -- Current landing pixel to crater center in real space
!   * xpi, ypi -- Current landing pixel to crater center in pixel space
!   * ebh      -- ejecta thickness
!   * rm       -- radius of melt zone
!   
!   Output
!   * surf     -- Output composition on a grid space 
! 
! NOTES
!   In future, a segment may contain multicomponent, and as a result the advanced analysis is needed. 
!
!***

!**********************************************************************************************************************************
!
!  Unit Name   : regolith_streamtube
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Determine a history of a stream tube
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  'eradc' is the center of the ejection radius. 'eradi' is inner, 'erado' is outer. 'cnt' is counting number. 'cmax', 'ri', and 'rip1' concern streamtube geometry.
!
!**********************************************************************************************************************************
subroutine regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,rm,vsq,volm)
   use module_globals 
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   integer(I4B),intent(in) :: ejtble
   type(ejbtype),dimension(:),intent(in)   :: ejb
   real(DP),intent(in)          :: xp,yp,lrad,ebh
   integer(I4B),intent(in)      :: xpi,ypi
   real(DP),intent(in)          :: rm, vsq
   real(DP),intent(inout)       :: volm

   ! Traversing a linked list 
   real(DP),parameter :: a = 0.936457 ! Fitting parameters for the relation between height difference and a radial position of a stream tube
   real(DP),parameter :: b = 1.12368
   !real(DP),parameter :: dz = 20.0
   real(DP)     :: frac,logtablerad,loglrad,logdelta,outeredge,inneredge
   real(DP)     :: deltar
   real(DP)     :: erado,eradi,xl,yl,eradc
   real(DP)     :: theta_eradi,length,vhead
   real(DP)     :: dy,ry,zo,zi
   real(DP),dimension(2) :: y   
   real(DP),dimension(:),allocatable :: xints,yints
   real(DP)     :: xi,yj,phi,cmax,tmpx,tmpy,rints
   integer(I4B) :: i,j,k,toti,totj,toty,cnt,xstpi,ystpi
   real(DP)     :: vtot,vseg,ri,rip1,xc,yc,thetast
   real(DP)     :: vst,vbody,rbody,vmare,totmare,totseb,tots
   real(DP)     :: meltinejecta, totvol, factor, agefactor
   type(regodatatype) :: newlayer
   real(SP),dimension(:),allocatable :: distvol

   ! Constrain the tangital tube's volume with CTEM result
   real(DP)     :: k1,k2,k3,k4,c1,c2

   ! check out if a tube exits a cell
   real(DP) :: dri,drip1
   integer(I4B) :: xril,xrir,yrib,yrit,xrip1l,xrip1r,yrip1b,yrip1t

   ! NaN and infinity number debugging 
   real(DP) :: x,xc1,xc2,yc1,yc2,xstpi1,xstpi2,ystpi1,ystpi2

   ! Monte Carlo method for two layers system
   !real(DP) :: compnohead,compmc,mceb
   ! Melt vs glass variables
   real(DP) :: meltfrac, melt, volm1, volv1, vol
   real(DP) :: xm, sinvints, cosvints, rvints, rmints, xvints
   real(DP) :: thetaq, depthb, q1, q2, q3, xmints

   ! Shock fragmentation variables
   real(DP) :: xsfints, rsh

   ! Melt vs age estimation
   real(SP),dimension(MAXAGEBINS) :: age_collector
   integer(I2B)               :: n_age

   ! Executalbe code

   ! ****** Interpolate radial distance, erad, for a given pixel *******
   ! outeredge = crater%frad + domain%ejbres * (EJBTABSIZE - 0.5_DP)
   ! inneredge = crater%frad + 0.5_DP * domain%ejbres
   ! k = max(min(1 + int((lrad - inneredge) / (outeredge - inneredge) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
   ! loglrad = log(lrad)
   ! logtablerad = ejb(k)%lrad 

   !from ejecta_interpolate
   inneredge = crater%ejrad 
   outeredge = crater%ejrad * exp(domain%ejbres * EJBTABSIZE)
   k = max(min(1 + int((log(lrad) - log(inneredge)) / (log(outeredge) - log(inneredge)) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
   loglrad = log(lrad)
   logtablerad = ejb(k)%lrad

   !from stable-1.4
   ! inneredge = crater%rad 
   ! outeredge = crater%rad * exp(domain%ejbres * EJBTABSIZE)
   ! k = max(min(1 + int((log(lrad) - log(inneredge)) / (log(outeredge) - log(inneredge)) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
   ! loglrad = log(lrad)
   ! logtablerad = ejb(k)%lrad 


   if (k==ejtble) then
      logdelta = logtablerad - ejb(k - 1)%lrad
      frac = (loglrad - ejb(k-1)%lrad) / logdelta
      eradc = exp(ejb(k-1)%erad) + ((exp(ejb(k)%erad) - exp(ejb(k-1)%erad)) * frac)
      !eradc = ejb(k-1)%erad + ((ejb(k)%erad - ejb(k-1)%erad) * frac)
   else
      logdelta = ejb(k + 1)%lrad - logtablerad 
      frac = (loglrad - logtablerad) / logdelta  
      eradc = exp(ejb(k)%erad) - ((exp(ejb(k)%erad) - exp(ejb(k+1)%erad)) * frac)
      !eradc = ejb(k)%erad - ((ejb(k)%erad - ejb(k+1)%erad) * frac)
   end if

   if (eradc<=0.0_DP) then
   write(*,*) k,ebh,crater%ejdis,lrad,exp(ejb(k-1)%lrad),exp(ejb(k)%lrad),eradc,exp(ejb(k-1)%erad),exp(ejb(k)%erad)
   stop
   end if
 
   ! ********* Calculate the height difference between two streamlines that define a stream tube with varying radial position *******
   xl = xp - crater%xl
   yl = yp - crater%yl
   phi = atan2(yl,xl) 
   toti = floor((eradc * abs(xl)/lrad) / user%pix)
   totj = floor((eradc * abs(yl)/lrad) / user%pix)
   ! *************************** Intersection points with grid lines ************************************************
   ! Allocate a space for total interection points with horizontal and vertical grid lines
   ! Finding intersection points starts from the the emerging point to the origin, so once we find all the intersection
   ! points, we can extract each segments intersected with horizontal grid line or vertical line and estimate the volume 
   ! of a segment of a stream tube by a scaled relationsip between the position of a segment and the width of a stream tube 
   ! at the position of a segment. 

   allocate(xints(abs(toti)+abs(totj) + 2))
   allocate(yints(abs(toti)+abs(totj) + 2))

   cnt = 1
   xints(cnt) = eradc * xl/lrad  
   yints(cnt) = eradc * yl/lrad 

   ! 1) Intersect with vertical grid line: start to calculate the intersection point from the larger absolute value of i or j
   !    , and the quadrant of an intersection point will be recovered from the sign of yl or xl that is shifted to the coordinate
   !    centered by the center of a crater. Later, it will be shited back to the grid space by adding the location of 
   !    a crater in pixel space. 

   do j=totj,1,-1
      yj = sign(1.0_DP,yl) * (real(j) * user%pix)
      cnt = cnt + 1
      xints(cnt) = yj/tan(phi)
      yints(cnt) = yj 
   end do

   do i=toti,1,-1
      xi = sign(1.0_DP,xl) * (real(i) * user%pix)
      cnt = cnt + 1
      xints(cnt) = xi
      yints(cnt) = tan(phi) * xi
   end do
 
   cnt = cnt + 1
   xints(cnt) = 0._DP
   yints(cnt) = 0._DP
 
   ! ******* Sort those final intersection points by a difference of a radial distance of an intersection point *********
   ! ******************************* Compare the radial distance of an intersected point ******************************
   do i=2,cnt-1
      cmax = sqrt(xints(i)**2 + yints(i)**2)
      do j=i+1,cnt-1
         rints = sqrt(xints(j)**2 + yints(j)**2)
         if (rints > cmax) then 
            cmax = rints  
            tmpx = xints(i)
            xints(i) = xints(j)
            xints(j) = tmpx
 
            tmpy = yints(i)
            yints(i) = yints(j)
            yints(j) = tmpy
         end if
      end do
   end do

   ! Purpose 1: Calculate the size of a stream tube constrained by ejecta thickness
   ! Purpose 2: Once we have the size information of a stream tube, we can
   ! calculate the distal melt: the precursor of glass spherules within a
   ! stream tube. The result is contained in a linked list "newlayer".
   call regolith_melt_glass(user,crater,domain,ebh,rm,eradc,lrad,deltar,newlayer,xmints,volm)
   ! if (eradc>rm) then
   !    write(*,*) 'eradc > rm!'
   !    write(*,*) ebh, exp(ejb(k)%thick)
   !   stop
   ! end if 
   erado = eradc + deltar
   eradi = eradc - deltar
   age_collector(:) = 0.0_SP
   vol = 0.0_DP
   totmare = 0.0_DP
   tots    = 0.0_DP
   depthb = crater%imp / 2.0_DP
   meltinejecta = 0.0_DP
   totvol = 0.0_DP
   allocate(distvol(1+domain%rcnum))
   distvol(:) = 0.0_SP

   ! if (eradc<=user%testimp) then
   !    write(*,*) lrad/crater%frad, user%testimp, crater%frad, rm, deltar, eradc, eradi, erado, ebh, newlayer%meltfrac
   !    stop
   ! end if

   call regolith_shock_damage_zone(crater,rm,eradi,depthb,xsfints)

   if (eradc <= user%pix) then

      xc = (xints(2) + xints(1))/2.0_DP
      yc = (yints(2) + yints(1))/2.0_DP
      xstpi = crater%xlpx + nint(xc/user%pix) 
      ystpi = crater%ylpx + nint(yc/user%pix) 
      ri = sqrt(xints(2)**2 + yints(2)**2)
      rip1 = sqrt(xints(1)**2 + yints(1)**2)
      vseg = regolith_streamtube_volume_func(eradi,0.0_DP,eradi,deltar)
      newlayer%thickness = vseg/(user%pix**2)
      call util_periodic(xstpi,ystpi,user%gridsize)
      call regolith_subpixel_streamtube(user,surf(xstpi,ystpi),deltar,ri,rip1,eradi,newlayer,vmare,totseb,&
           age_collector,xmints,xsfints,vol,meltinejecta,totvol,distvol)

      newlayer%age(:) = newlayer%age(:) + age_collector(:)
 
      totmare = vmare
      tots = totseb
      newlayer%thickness = ebh
      newlayer%comp      = min(totmare/tots, 1.0_DP)
      !newlayer%age(:)    = newlayer%age(:) * min( (ebh * user%pix**2) / tots, 1.0_DP)

   else
      rbody = sqrt(xints(2)**2 + yints(2)**2)
      if (rbody<eradi) then
         xc = (xints(2) + eradi*xl/lrad)/2.0
         yc = (yints(2) + eradi*yl/lrad)/2.0
         xstpi = crater%xlpx + nint(xc/user%pix)
         ystpi = crater%ylpx + nint(yc/user%pix)
         vseg = regolith_streamtube_volume_func(eradi,rbody,eradi,deltar)
         newlayer%thickness = vseg/(user%pix**2)
         call util_periodic(xstpi,ystpi,user%gridsize)
         call regolith_traverse_streamtube(user,surf(xstpi,ystpi),deltar,rbody,eradi,eradi,erado,newlayer,vmare,&
            totseb,age_collector,xmints,xsfints,depthb,meltinejecta,totvol,distvol)
         totmare = totmare + vmare
         tots = tots + totseb
      end if           
  
      do i=cnt,3,-1
         xc = (xints(i) + xints(i-1))/2.0_DP
         yc = (yints(i) + xints(i-1))/2.0_DP
         xstpi = crater%xlpx + nint(xc/user%pix)
         ystpi = crater%ylpx + nint(yc/user%pix) 
         ri = sqrt(xints(i)**2 + yints(i)**2)
         rip1 = sqrt(xints(i-1)**2 + yints(i-1)**2)
         if (abs(ri-rip1)>user%pix/1000000.0) then
            vseg = regolith_streamtube_volume_func(eradi,ri,rip1,deltar)
            newlayer%thickness = vseg/(user%pix**2)
            call util_periodic(xstpi,ystpi,user%gridsize)
            call regolith_traverse_streamtube(user,surf(xstpi,ystpi),deltar,ri,rip1,eradi,erado,newlayer,vmare,&
               totseb,age_collector,xmints,xsfints,depthb,meltinejecta,totvol,distvol)
            totmare = totmare + vmare
            tots = tots + totseb 
         end if
      end do

      xstpi = crater%xlpx + nint(eradc*xl/lrad/user%pix)
      ystpi = crater%ylpx + nint(eradc*yl/lrad/user%pix)
      call util_periodic(xstpi,ystpi,user%gridsize)
      call regolith_streamtube_head(user,surf(xstpi,ystpi),deltar,totmare,tots,age_collector,meltinejecta,totvol,distvol)

      newlayer%thickness = ebh
      newlayer%comp      = min(totmare/tots, 1.0_DP)
      newlayer%age(:)    = newlayer%age(:) + age_collector(:)
      !newlayer%age(:)    = newlayer%age(:) * min( (ebh * user%pix**2) / tots, 1.0_DP)
      ! if (newlayer%meltfrac > 1.0_DP) then
      !    write(*,*) "Melt fraction >1! (Traverse)", xpi,ypi,crater%timestamp,crater%fcrat,crater%xlpx,crater%ylpx,&
      !     newlayer%meltvolume, newlayer%totvolume, newlayer%ejm, newlayer%ejmf, totvol
      ! end if
   end if

  !Apply a correction factor to ensure conservation of volume

  factor = (newlayer%totvolume-newlayer%ejm) / totvol
  meltinejecta = meltinejecta * factor
  distvol(:) = distvol(:) * factor
  !totvol = newlayer%totvolume - meltinejecta
  if (newlayer%ejm > newlayer%totvolume) then !entire pixel is ejected melt
      newlayer%ejm = newlayer%totvolume
      newlayer%meltvolume = newlayer%ejm
   else
      if (meltinejecta + newlayer%ejm > newlayer%totvolume) then !entire pixel is melt, but not all of it is ejected
         meltinejecta = newlayer%totvolume - newlayer%ejm
      end if
      newlayer%meltvolume = meltinejecta + newlayer%ejm
      if (newlayer%meltvolume > newlayer%totvolume) then !edge case caused by floating point math could result in melt fraction slightly higher than 1
         factor = newlayer%totvolume / newlayer%meltvolume
         newlayer%meltvolume = newlayer%totvolume
         distvol(:) = distvol(:) * factor
         newlayer%age(:) = newlayer%age(:) * factor
      end if
      newlayer%distvol(:) = newlayer%distvol(:) + distvol(:)

      newlayer%distvol(1+domain%rcnum) = newlayer%meltvolume - sum(newlayer%distvol(1:domain%rcnum))
      if (newlayer%distvol(1+domain%rcnum) < 0.0) then !pixel consists entirely of QMC melt
         newlayer%distvol(1+domain%rcnum) = 0.0_SP
         newlayer%age(:) = 0.0_SP
      end if
      if (sum(newlayer%distvol) > newlayer%totvolume) then
         factor = newlayer%totvolume / sum(newlayer%distvol)
         newlayer%distvol(:) = newlayer%distvol(:) * factor
         newlayer%age(:) = newlayer%age(:) * factor

      end if
      newlayer%meltvolume = sum(newlayer%distvol)

   end if

   !conserve volume in the age array
   if (sum(newlayer%age(:)) > 0.0) then
      agefactor = newlayer%distvol(1+domain%rcnum) / sum(newlayer%age(:))
      newlayer%age(:) = newlayer%age(:) * agefactor
   else
      newlayer%age(:) = 0.0_SP
   end if

   if (.not. allocated(newlayer%regotemp)) then
      allocate(newlayer%regotemp(1,1))
      newlayer%regotemp(1,1) = 0.0_SP
   end if
   if (.not. allocated(newlayer%regotime)) then
      allocate(newlayer%regotime(1,1))
      newlayer%regotime(1,1) = 0.0_SP
   end if

  call util_push_array(surf(xpi,ypi)%regolayer,newlayer)

  deallocate(xints,yints,distvol)

  return
end subroutine regolith_streamtube
