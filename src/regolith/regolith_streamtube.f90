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
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_streamtube(user,surf,crater,domain,ejb,ejtble,xp,yp,xpi,ypi,lrad,ebh,rm)
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
   real(DP),intent(in)          :: rm

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
   type(regodatatype) :: newlayer, sumlayer

   ! Constrain the tangital tube's volume with CTEM result
   real(DP)     :: k1,k2,k3,k4,c1,c2

   ! check out if a tube exits a cell
   real(DP) :: dri,drip1
   integer(I4B) :: xril,xrir,yrib,yrit,xrip1l,xrip1r,yrip1b,yrip1t

   ! NaN and infinity number debugging 
   real(DP) :: x,xc1,xc2,yc1,yc2,xstpi1,xstpi2,ystpi1,ystpi2

   ! Monte Carlo method for two layers system
   !real(DP) :: compnohead,compmc,mceb

   ! Fresh melt's calculation
   real(DP) :: meltfrac, melt, volm1, volv1
   real(DP) :: xm, sinvints, cosvints, rvints, rmints, xvints
   real(DP) :: cosq1, cosq2, thetaq, depthb, q1, q2, q3, xmints   

   ! Executalbe code

   ! ****** Interpolate radial distance, erad, for a given pixel *******
   outeredge = crater%frad + domain%ejbres * (EJBTABSIZE - 0.5_DP)
   inneredge = crater%frad + 0.5_DP * domain%ejbres
   k = max(min(1 + int((lrad - inneredge) / (outeredge - inneredge) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
   loglrad = log(lrad)
   logtablerad = ejb(k)%lrad 

   if (k==ejtble) then
      logdelta = logtablerad - ejb(k - 1)%lrad
      frac = (loglrad - ejb(k-1)%lrad) / logdelta
      eradc = exp(ejb(k-1)%erad) + ((exp(ejb(k)%erad) - exp(ejb(k-1)%erad)) * frac)
   !   frac = (loglrad - logtablerad) / logdelta
   !   eradc = ejb(k)%erad - ((ejb(k)%erad - LOGVSMALL) * frac)
   else
      logdelta = ejb(k + 1)%lrad - logtablerad 
      frac = (loglrad - logtablerad) / logdelta  
      eradc = exp(ejb(k)%erad) - ((exp(ejb(k)%erad) - exp(ejb(k+1)%erad)) * frac)
   end if

   if (eradc<=0.0_DP) then
   write(*,*) k,ebh,crater%ejdis,lrad,exp(ejb(k-1)%lrad),exp(ejb(k)%lrad),eradc,ejb(k-1)%erad,ejb(k)%erad
   stop
   end if

   ! ********* Calculate the height difference between two streamlines that define a stream tube with varying radial position *******
   xl = xp - crater%xl
   yl = yp - crater%yl
   phi = atan(yl/xl) 
   !toti = ceiling((eradc * abs(xl)/lrad - 0.5_DP * user%pix) / user%pix) ! Refer to the definition of CTEM in ejecta_emplace.f90
   !totj = ceiling((eradc * abs(yl)/lrad - 0.5_DP * user%pix) / user%pix) !
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

   !write(*,*) toti,totj
   ! 1) Intersect with vertical grid line: start to calculate the intersection point from the larger absolute value of i or j
   !    , and the quadrant of an intersection point will be recovered from the sign of yl or xl that is shifted to the coordinate
   !    centered by the center of a crater. Later, it will be shited back to the grid space by adding the location of 
   !    a crater in pixel space. 

   do j=totj,1,-1
      !yj = sign(1.0_DP,yl) * (real(j-1) * user%pix + 0.5_DP * user%pix) 
      yj = sign(1.0_DP,yl) * (real(j) * user%pix)
      cnt = cnt + 1
      xints(cnt) = yj/tan(phi)
      yints(cnt) = yj 
      !if (abs(phi) < VSMALL .or. (abs(phi) >= PI/2.0_DP - VSMALL .and. abs(phi) <= PI/2.0_DP + VSMALL) ) then
      !if (totj == 0) then
      !   write(*,*) 'y',phi,cnt,yj,xints(cnt),yints(cnt),sqrt(xints(cnt)**2 + yints(cnt)**2)
      !end if
   end do

   do i=toti,1,-1
      !xi = sign(1.0_DP,xl) * (real(i-1) * user%pix + 0.5_DP * user%pix)
      xi = sign(1.0_DP,xl) * (real(i) * user%pix)
      cnt = cnt + 1
      xints(cnt) = xi
      yints(cnt) = tan(phi) * xi
      !if (abs(phi) < VSMALL .or. (abs(phi) >= PI/2.0_DP - VSMALL .and. abs(phi) <= PI/2.0_DP + VSMALL) ) then
      !   write(*,*) 'x',phi,cnt,xi,xints(cnt),yints(cnt),sqrt(xints(cnt)**2 + yints(cnt)**2)
      !end if
   end do
 
   cnt = cnt + 1
   xints(cnt) = 0._DP
   yints(cnt) = 0._DP

   !if (lrad/1000.0<=1.0 .and. lrad/1000.0>=0.98) then
   !   do i=cnt,1,-1
   !      write(*,*) lrad/1000.0,xints(i),yints(i),sqrt(xints(i)**2 + yints(i)**2)
   !   end do
   !end if
 
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

   ! ******************************** Size of a tangential stream tube **********************************************
   ! Determine the radius of a tangential-shaped circule stream tube at the head of the stream tube
   ! The ejecta blanket thickness data from CTEM will be used to constrain our tube, then two calibrated 
   ! streamlines will be given for estimation of thickness of a layer inside a stream tube
   ! A cubic function is given after considering the constraint of ejecta blanket thickness by CTEM, so the cubic 
   ! function's solution is based on Nmuerical Recipe: Fortran 77, p.178-200, which is believed in a optimized way to do so!
   ! The total volume of a stream tube is given:
   ! Vst = (0.25_DP * PI * deltar**2 * a**2 * eradi / b *(tan(b)-b)) + sqrt(2.0_DP)/2.0_DP*PI*deltar**3
   
   k1 = PI * (0.5 * a)**2 * (tan(b) - b)
   k2 = ebh * (user%pix)**2
   k3 = eradc/b 
   k4 = sqrt(3.0_DP)/2.0_DP * PI - k1/b
   c1 = k1 * k3 /k4
   c2 = k2/k4
   deltar = regolith_cubic_func(c1,c2)
   erado = eradc + deltar
   eradi = eradc - deltar
   !if (deltar <= VSMALL) write(*,*) crater%frad,lrad/crater%frad,eradc,deltar,ebh,k1,k2,k3,k4,c1,c2
   
   ! Comparison with melt zone's radius
   !write(*,*) lrad / crater%frad, eradi, erado, deltar


   ! ******************************* Start to estimate STREAM TUBE'S volume in layering systerm ***************************
   ! Purpose: How much layer material are contained in a stream tube? 
   ! Intro: There are two volume approximation with regarding to discretized stream tubes. First, the subpixel approximation 
   ! is important to distal landing ejecta, and the advantage of it is that you can use Maxwell Z model equation to calculate
   ! the total volume of a stream tube if it is contained inside the whole layer. For the small craters or subpixel craters,
   ! this method will be used frequently. 

   vst = (0.25_DP * PI * deltar**2 * a**2 * eradi / b *(tan(b)-b)) + sqrt(2.0_DP)/2.0_DP*PI*deltar**3

   ! Fresh melt's calculation
   depthb             = crater%imp / 2.0
   rints              = sqrt(rm**2 - (depthb)**2)
   cosvints           = eradi / (crater%imp + eradi)
   sinvints           = sqrt(1.0 - cosvints**2)
   xvints             = eradi * ( 1.0 - cosvints) * sinvints
   volv1              = ( 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (tan(b/eradi*(xvints))-b/eradi*(xvints)) )
   if (eradi <= rints) then
      volm1             = vst - volv1
      melt              = volm1
      newlayer%meltfrac = 1.0
   else if (eradi > rints) then
           q1           =   1.0 / (1.0 + 2.0 * depthb / eradi)
           q2           =  -1.0 - 1.0/q1
           q3           = ( 1.0 + (depthb**2 - rm**2)/eradi**2 ) * q1
           cosq1        = 0.5 * q1 * (-1.0 * q2 + sqrt(q2**2 - 4.0 * q3/q1))
           cosq2        = 0.5 * q1 * (-1.0 * q2 - sqrt(q2**2 - 4.0 * q3/q1))
           thetaq       = acos( min(abs(cosq1),abs(cosq2)) )
           xmints       = eradi * (1.0 - cos(thetaq)) * sin(thetaq)
           volm1        = ( 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (tan(b/eradi*(xmints))-b/eradi*(xmints)) )
           melt         = volm1 - volv1
           newlayer%meltfrac = melt/(vst-volv1)
   end if

   if (eradc <= user%pix) then

      xc = (xints(2) + xints(1))/2.0_DP
      yc = (yints(2) + yints(1))/2.0_DP
      xstpi = crater%xlpx + nint(xc/user%pix) 
      ystpi = crater%ylpx + nint(yc/user%pix) 
      ri = sqrt(xints(2)**2 + yints(2)**2)
      rip1 = sqrt(xints(1)**2 + yints(1)**2)
      vseg = 0.25_DP * PI * deltar**2 * a**2 * eradi / b * abs(tan(b) - b) 
      call util_periodic(xstpi,ystpi,user%gridsize)
      call regolith_subpixel_streamtube(user,surf(xstpi,ystpi),deltar,ri,rip1,eradi,vseg,newlayer,rm)
      newlayer%comp = newlayer%comp / newlayer%thickness
      newlayer%thickness = ebh

   else
      ! Define a regodata typed layer that is able to sum up all components during travesing a stream tube
      sumlayer%thickness = 0._DP
      sumlayer%comp      = 0._DP
      sumlayer%meltfrac  = 0._DP 

      rbody = sqrt(xints(2)**2 + yints(2)**2)
      if (rbody<eradi) then
         xc = (xints(2) + eradi*xl/lrad)/2.0
         yc = (yints(2) + eradi*yl/lrad)/2.0
         xstpi = crater%xlpx + nint(xc/user%pix)
         ystpi = crater%ylpx + nint(yc/user%pix)
         vseg = 0.25 * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b) - tan(b/eradi * rbody)) - &
               abs(b - b/eradi * rbody))
         call util_periodic(xstpi,ystpi,user%gridsize)
         call regolith_traverse_streamtube(user,surf(xstpi,ystpi),deltar,rbody,eradi,eradi,erado,vseg,newlayer,rm)
         sumlayer%comp      = sumlayer%comp      + newlayer%comp
         sumlayer%thickness = sumlayer%thickness + newlayer%thickness
      end if           
  
      do i=cnt,3,-1
         xc = (xints(i) + xints(i-1))/2.0_DP
         yc = (yints(i) + xints(i-1))/2.0_DP
         xstpi = crater%xlpx + nint(xc/user%pix)
         ystpi = crater%ylpx + nint(yc/user%pix) 
         ri = sqrt(xints(i)**2 + yints(i)**2)
         rip1 = sqrt(xints(i-1)**2 + yints(i-1)**2)
         if (abs(ri-rip1)>user%pix/1000000.0) then
            vseg = 0.25_DP * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b/eradi * rip1) - tan(b/eradi * ri)) &
                   - abs(b/eradi * rip1 - b/eradi * ri)) 
            call util_periodic(xstpi,ystpi,user%gridsize)
            call regolith_traverse_streamtube(user,surf(xstpi,ystpi),deltar,ri,rip1,eradi,erado,vseg,newlayer,rm)
            sumlayer%comp      = sumlayer%comp      + newlayer%comp
            sumlayer%thickness = sumlayer%thickness + newlayer%thickness
         end if
      end do

      newlayer%thickness = 0._DP
      newlayer%comp      = 0._DP
      xstpi = crater%xlpx + nint(eradc*xl/lrad/user%pix)
      ystpi = crater%ylpx + nint(eradc*yl/lrad/user%pix)
      call util_periodic(xstpi,ystpi,user%gridsize)
      call regolith_streamtube_head(user,surf(xstpi,ystpi),deltar,newlayer,eradi,rm)
      sumlayer%comp      = sumlayer%comp      + newlayer%comp
      sumlayer%thickness = sumlayer%thickness + newlayer%thickness
      newlayer%comp      = sumlayer%comp / sumlayer%thickness
      newlayer%thickness = ebh

  end if

  call util_push(surf(xpi,ypi)%regolayer,newlayer)
  deallocate(xints,yints)

  return
end subroutine regolith_streamtube
