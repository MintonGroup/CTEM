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
   use module_util
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
