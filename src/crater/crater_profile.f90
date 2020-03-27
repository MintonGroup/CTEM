!**********************************************************************************************************************************
!
!  Unit Name   : crater_profile
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Function that defines the basic profile of a crater, not including the ejecta blanket
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
function crater_profile(user,crater,r) result(h)
   use module_globals
   use module_crater, EXCEPT_THIS_ONE => crater_profile
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: r
   real(DP) :: h

   ! Internal variables
   real(DP) :: flrad,c0,c1,c2,c3

   ! Executable code
   flrad = 0.5_DP * crater%floordiam / crater%frad 
   
   ! Use polynomial crater profile similar to that of Fassett et al. (2014), but the parameters are set by the crater dimensions
   c1 = (-crater%floordepth - crater%rimheight) / (flrad + flrad**2 / 3._DP - flrad**3 / 6._DP - 7._DP / 6._DP)
   c0 = crater%rimheight - (7._DP / 6._DP) * c1
   c2 = c1 / 3._DP
   c3 = -c2 / 2._DP

   if (r < flrad) then
      h = -crater%floordepth 
   if (r > crater%frad)
      h = crater%rimheight * ((crater%frad / lrad)**RIMDROP)
   else
      h = c0 + c1 * r + c2 * r**2 + c3 * r**3 
   end if


   return
end subroutine crater_profile

