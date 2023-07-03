!****f* regolith/regolith_melt_glass
! Name
!   regolith_melt_glass -- Calculate melt glass within a stream tube during excavation stage.
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_regolith
!   
!   call regolith_melt_glass(user,age,ebh,deltar,newlayer)
!
! DESCRIPTION
!
!   The modeling of production of glass spherules is based on the terrestrial
!   tektites and microtektites record. A direct cutoff (no glass spherules
!   passing a specific landing distance) is applied. Future work can be done to
!   interpret the efficiency of glass spherule production from their trajectory
!   from iSale and check interaction between hot ejecta, vapor, and cold ejecta.
!
! ARGUMENTS
!   Input
!   * user     -- The user-defined variables from the input file
!   * age      -- The age of modeled glass spherule which is equivalent to the
!                 formtion of its parent impact crater.
!   * ebh      -- The thickness of ejecta blanket
!
!   Output
!   * deltar   -- The size or cross section's radius of a stream tube
!   * newlayer -- The linked list that contains ejecta's properites.
!   
! NOTES
!   In future, a segment may contain multicomponent, and as a result the advanced analysis is needed. 
!
!***

!**********************************************************************************************************************************
!
!  Unit Name   : regolith_melt_glass
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Determine an abundance of glass spherules within a stream tube
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : 
!      newlayer 
!      deltar
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_melt_glass(user,crater,domain,age,age_resolution,ebh,rm,eradc,lrad,deltar,newlayer,xmints,melt)
   use module_globals 
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_melt_glass
   implicit none

   ! Arguments
   type(usertype),intent(in)        :: user
   type(cratertype),intent(in)      :: crater
   type(domaintype),intent(in)      :: domain
   real(DP),intent(in)              :: age
   real(DP),intent(in)              :: age_resolution
   real(DP),intent(in)              :: ebh
   real(DP),intent(in)              :: rm
   real(DP),intent(in)              :: eradc
   real(DP),intent(in)              :: lrad
   real(DP),intent(out)             :: deltar
   type(regodatatype),intent(out)   :: newlayer
   real(DP),intent(out)             :: xmints
   real(DP),intent(out)             :: melt

   ! Internal variables
   ! Stream tube parameters  
   real(DP),parameter :: a = 0.936457 ! Fitting parameters for the relation between height difference and a radial position of a stream tube
   real(DP),parameter :: b = 1.12368

   ! Constrain the tangital tube's volume with CTEM result
   real(DP)     :: k1,k2,k3,k4,c1,c2
   
   ! Calculate vapor and melt zone intersection point with stream tubes
   real(DP)     :: vst, erado, eradi
   real(DP)     :: cosvints, sinvints, xvints, rints
   real(DP)     :: volv1, volm1, depthb
   real(DP)     :: q1, q2, q3 
   real(DP)     :: thetaq
   integer(I2B) :: n_age
   ! Estimate melt droplet size
   real(DP),parameter :: b_exponent = -0.97
   real(DP)           :: cvpgsqr, p1, p2, p3, p4, p5
   real(DP)           :: dm
   real(DP)           :: cosq1, cosq2
   

   ! Executalbe code
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

   ! ******************************* Start to estimate STREAM TUBE'S volume in layering systerm ***************************
   ! Purpose: How much layer material are contained in a stream tube? 
   ! Intro: There are two volume approximation with regarding to discretized stream tubes. First, the subpixel approximation 
   ! is important to distal landing ejecta, and the advantage of it is that you can use Maxwell Z model equation to calculate
   ! the total volume of a stream tube if it is contained inside the whole layer. For the small craters or subpixel craters,
   ! this method will be used frequently. 

   vst = (0.25_DP * PI * deltar**2 * a**2 * eradi / b *(tan(b)-b)) + sqrt(2.0_DP)/2.0_DP*PI*deltar**3
   
   newlayer%thickness = ebh    ! default value: stream tube's volume = paraboloid shell's volume
   newlayer%comp      = 0.0_DP
   newlayer%meltvolume = 0.0_DP
   newlayer%totvolume = newlayer%thickness * user%pix * user%pix
   newlayer%ejm       = 0.0_DP
   rints              = sqrt(rm**2 - (crater%imp/2.0)**2)
   cosvints           = min(max(eradi / (crater%imp + eradi), -1.0_DP), 1.0_DP)
   sinvints           = sqrt(1.0 - cosvints**2)
   xvints             = eradi * ( 1.0 - cosvints) * sinvints
   volv1 = regolith_streamtube_volume_func(eradi,0.0_DP,xvints,deltar)
   melt               = 0.0_DP
   newlayer%age(:)    = 0.0_SP
   if (eradi <= rints) then
      volm1    = vst - volv1
      melt     = volm1
      newlayer%meltvolume = melt
      !newlayer%totvolume = volm1
      newlayer%ejm = melt
      xmints   = rints 
   else if (eradi > rints) then
           depthb = crater%imp / 2.0
         !   q1     =   1.0 / (1.0 + 2.0 * depthb / eradi)
         !   q2     =  -1.0 - q1
         !   q3     = ( 1.0 + (depthb**2 - rm**2)/eradi**2 ) * q1 
         !   thetaq = acos( -0.5 * q2 - 0.5 * sqrt(q2**2 - 4.0 * q3) )
         !   xmints = eradi * (1.0 - cos(thetaq)) * sin(thetaq)
         !   volm1 = regolith_streamtube_volume_func(eradi,0.0_DP,xmints,deltar)
         !   melt   = volm1 - volv1
         !   newlayer%meltfrac = melt/vst

           !the following is from the old regolith_streamtube.f90:

           q1           =   1.0 / (1.0 + 2.0 * depthb / eradi)
           q2           =  -1.0 - 1.0/q1
           q3           = ( 1.0 + (depthb**2 - rm**2)/eradi**2 ) * q1
           cosq1        = 0.5 * q1 * (-1.0 * q2 + sqrt(q2**2 - 4.0 * q3/q1))
           cosq2        = 0.5 * q1 * (-1.0 * q2 - sqrt(q2**2 - 4.0 * q3/q1))
           thetaq       = acos( min(abs(cosq1),abs(cosq2)) )
           xmints       = eradi * (1.0 - cos(thetaq)) * sin(thetaq)
           volm1        = regolith_streamtube_volume_func(eradi,0.0_DP,xmints,deltar)
           melt         = volm1 - volv1
           newlayer%meltvolume = melt
           newlayer%totvolume = newlayer%thickness * user%pix * user%pix
           newlayer%ejm = melt
           
   end if

   allocate(newlayer%distvol((1+domain%rcnum)))
   newlayer%distvol(:) = 0.0_SP
   if(domain%currentqmc) then
      newlayer%distvol(domain%nqmc) = newlayer%meltvolume
   else
      newlayer%distvol(1+domain%rcnum) = newlayer%meltvolume
      newlayer%age(domain%age_counter) = newlayer%meltvolume
   end if 

   ! n_age = max(ceiling(age / age_resolution), 1)
   ! if (lrad >= RAD_GP * crater%rad) then
   !    newlayer%age(n_age) = melt / (user%pix * user%pix) 
   ! else 
   !    newlayer%age(n_age) = 0.0_SP   
   ! end if

   return
end subroutine regolith_melt_glass
