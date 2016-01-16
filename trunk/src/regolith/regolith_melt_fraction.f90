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
   
   maxerad = max(erad1,erad2)
   minerad = min(erad1,erad2)
   rvapor  = 0.5 * dimp 
   rints = sqrt(rmelt**2 - depthb**2)

   if (maxerad<=rints) then 
           meltfrac   = 1.0
   else if (minerad >= rints) then 
           volv1  = regolith_melt_func(rvapor,depthb,erad1)
           volv2  = regolith_melt_func(rvapor,depthb,erad2)           
           volm1  = regolith_melt_func(rmelt,depthb,erad1) - volv1
           volm2  = regolith_melt_func(rmelt,depthb,erad2) - volv2
           volm   = abs(volm2 - volm1)
           volst  = ( PI/6.0_DP * abs(erad1**3 - erad2**3) - abs(volv1 - volv2) ) 
           meltfrac = volm/volst
   end if

   return
end subroutine regolith_melt_fraction
