!**********************************************************************************************************************************
!
!  Unit Name   : regolith_melt_zone_superdomain
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Estimate the melt zone size of a superdomain crater
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
subroutine regolith_melt_zone_superdomain(user,crater,domain,rm,depthb)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_melt_zone_superdomain
   implicit none

   ! Arguments
   type(usertype),intent(in)             :: user
   type(cratertype),intent(inout)        :: crater
   type(domaintype),intent(in)           :: domain
   real(DP),intent(out)                  :: rm, depthb

   ! Internal variables
   real(DP)                     :: cvpg
   real(DP)                     :: vimp, sinimp, rimp, dimp
   real(DP)                     :: pvelv, pmass
   real(DP)                     :: c1, c2, pitwo, pithree, pifour, pivolg

      
   crater%grad   = crater%rad
   cvpg          = sqrt(2.0_DP) / 0.85_DP * ( user%mu_r / (user%mu_r + 1) )
   vimp          = domain%rmsvel
   crater%impvel = vimp
   crater%sinimpang = sqrt(2.0_DP) / 2.0_DP
   pvelv         = vimp * crater%sinimpang 
   pmass         = 4 * THIRD * PI * user%prho 
   c1            = 1.0_DP + 0.5_DP * user%mu_b
   c2            = (-3 * user%mu_r) / (2.0 + user%mu_r)
   pitwo         = (user%gaccel / pvelv**2)
   pifour        = user%trho_r / user%prho
   rimp          = ( PI/(3.0 * user%kv_r) * (crater%grad)**3 * (pitwo**(-1.0 * c2)) * &
                    (pifour**(c2/3.0)) * user%trho_r / pmass )**(1.0 / (3.0 + c2))
   dimp          = rimp * 2.0_DP 
   crater%imp    = dimp
   call regolith_melt_zone(user,crater,dimp,vimp,rm,depthb)
   return
end subroutine regolith_melt_zone_superdomain
