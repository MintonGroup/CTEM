!****f* regolith/regolith_melt_fraction
! Name
!   regolith_melt_fraction -- Calculate melt fraction of a stream tube 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_regolith
!   
!   call regolith_melt_fraction(dimp,depthb,erad1,erad2,rmelt,meltfrac)
!
! DESCRIPTION
!    
!   This subroutine estimates a ratio of melt to the volume of a stream tube. One must be cautious about the method that
!   we use to get a melt fraction. We assume that melt zone within a transient crater is homogeneously distributed. It is
!   the same in any direction from an impact site. We took this advantage and *revolution* a stream tube superposed with 
!   a vapor and melt zone. For example, two streamlines must overlap with a vapor and melt zone. To obtain melt proportion,
!   we need to know a vapor proportion and the total volume that two streamlines revolution. We revolution the function that
!   determines the one of two streamlines and the vapor zone and obtain the revolution volume. Then, we do the same thing 
!   for the second streamline. At the end, the substract volume between these two revolution volume is for vapor proportion 
!   inside these two streamlines zone. In more detail, there are three situation about intersection between two streamlines 
!   and vapor/melt zone: 1) two streamlines are completely inside melt zone. 2) two streamlines are completely outside melt 
!   zone. 3) two streamlines are inbetween outside and inside of melt zone. 
! 
! ARGUMENTS
!   Input
!   * dimp     -- The diameter of an impactor
!   * depthb   -- The burial depth of a melt zone's center relaitve to the surface
!   * erad1    -- A first streamline's emerging place from an impact site
!   * erad2    -- A second streamline's emerging place from an impact site
!   * rmelt    -- The radius of a melt zone
!
!   Output
!   * meltfrac -- Output a final calculation of melt fraction 
! 
! Notes
!   We do not consider streamlines inside a vapor zone.  
!
!***

!**********************************************************************************************************************************
!
!  Unit Name   : regolith_melt_fraction
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculate melt fraction in two adjacent streamlines 
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
subroutine regolith_melt_fraction(dimp,depthb,erad1,erad2,rmelt,meltfrac)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_melt_fraction
   implicit none

   ! Arguments
   real(DP),intent(in)       :: dimp,depthb,erad1,erad2,rmelt
   real(DP),intent(inout)    :: meltfrac

   ! Internal variables
   real(DP)                :: volm1,volm2,volm
   real(DP)                :: volst,volv1,volv2
   real(DP)                :: maxerad,minerad,rvapor,rints
   real(DP)                :: volv, volmi, volmo
   
   maxerad = max(erad1,erad2)
   minerad = min(erad1,erad2)
   rvapor  = 0.5 * dimp 
   rints = sqrt(rmelt**2 - depthb**2)

   if (maxerad<=rints) then 
           meltfrac   = 1.0
   else if (minerad >= rints) then 
           volv1  = regolith_melt_func(rvapor,depthb,minerad)
           volv2  = regolith_melt_func(rvapor,depthb,maxerad)           
           volm1  = regolith_melt_func(rmelt,depthb,minerad) - volv1
           volm2  = regolith_melt_func(rmelt,depthb,maxerad) - volv2
           volm   = abs(volm2 - volm1)
           volst  = ( PI/6.0_DP * abs(maxerad**3 - minerad**3) - abs(volv1 - volv2) ) 
           meltfrac = volm/volst
   else if (maxerad > rints .and. minerad < rints) then
           ! Vapor part inside a stream tube
           volv1  = regolith_melt_func(rvapor,depthb,minerad)
           volv2  = regolith_melt_func(rvapor,depthb,maxerad)
           volv   = abs(volv2 - volv1)
           ! Inside melt zone part
           volmi  = PI/6.0 * (rints**3 - minerad**3)
           volmo  = regolith_melt_func(rmelt,depthb,maxerad) - PI/6.0 * (rints)**3
           volm   = volmo + volmi - volv
           volst  = ( PI/6.0_DP * abs(maxerad**3 - minerad**3) - volv )
           meltfrac = volm/volst
   else
     write(*,*) 'regolith_melt_fraction: this is a bug!'
   end if

   return
end subroutine regolith_melt_fraction
