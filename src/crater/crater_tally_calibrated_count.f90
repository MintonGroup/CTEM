!**********************************************************************************************************************************
!
!  Unit Name   : crater_tally_calibrated_count
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Uses a human-calibrated counting model to determine the  countability of a given crater
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
subroutine crater_tally_calibrated_count(user,diameter,current_depth,original_depth,deviation_sigma,&
                                         countable,killable,p)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_tally_calibrated_count
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   real(DP),intent(in) :: diameter
   real(SP),intent(in) :: current_depth,original_depth,deviation_sigma
   logical,intent(out) :: countable,killable
   real(SP),intent(out) :: p

   ! Internal variables
   real(DP)               :: Rd,Rsig,pd,psig,a,b,c,d,Dtran,complex_correction

   Rd = current_depth/original_depth
   Rsig = log10(deviation_sigma/current_depth)

   ! Executable code
   select case(user%countingmodel)
   case("FASSETT")
      a = -0.0237520259849034_DP
      b = 0.208739267563691_DP
      c = 0.992603771065306_DP
      if (Rd < 0.081895370584072_DP) then
         pd = 0._DP
      else if (Rd > 0.9158503765541_DP) then
         pd = 1._DP
      else
         pd = a + b * Rd + c * Rd**2
      end if
      pd = min(max(pd, 0.0_DP), 1._DP)
 
      a = -0.0856409196596715_DP
      b = -1.16657321813299_DP
      psig = a + b * Rsig 
      psig = min(max(psig, 0.0_DP), 1._DP)

      p = real(min(pd,psig),kind = SP)
      p = min(max(p, 0.0_SP), 1._SP)
   case("HOWL")
         a = -0.48_DP
         b = -0.55_DP
         c = 8.0_DP 
         d = -4.5_DP
         p = (Rsig - a * log(Rd)) / b - c * (diameter / user%pix)**d - 0.5_DP
   end select

   select case (user%mat)
   case ("ROCK")
      Dtran = 2 * SIMCOMKS * user%gaccel**SIMCOMPS
   case ("ICE")
      Dtran = 2 * SIMCOMKI * user%gaccel**SIMCOMPI
   case default
      Dtran = 2 * SIMCOMKS * user%gaccel**SIMCOMPS ! Use rock as default just in case it dien't get set
   end select

   complex_correction = min(Dtran / diameter, 1.0_DP)

   if (p > 0.50_SP * complex_correction) then ! Kill the crater if its p-score is too low
      countable = .true.
      killable = .false.
   else
      countable = .false.
      killable = .true.
   end if

   ! TESTING
   if (user%testtally) then
      countable = .true.
      killable = .false.
   end if

   return
end subroutine crater_tally_calibrated_count

