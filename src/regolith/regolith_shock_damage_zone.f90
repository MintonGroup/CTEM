!****f* regolith/regolith_shock_damage_zone
! Name
!   regolith_shock_damage -- Calculate the radius of shock damage zone
!   and the intersection point between shock damage zone and the stream
!   tube
! SYNOPSIS
!   This uses
!   * module_globals
!   * module_regolith
!
!   call regolith_shock_damage_zone(user,crater,rm,eradi,depthb,xsfints)
!
! DESCRIPTION
!   
!   This subroutine calculates the radius of shock damage zone by a
!   given shock pressure, setting from globals module, and returns the
!   intersection point's x location between shock zone and a given
!   stream tube.
!  
! ARGUMENTS
!   Input
!   * user      -- The user-defined variables from the input file 
!   * crater    -- Crater dimension container
!   * rm        -- Radius of melt zone
!   * eradi     -- The inner radial distance of a stream tube
!   * depthb    -- The burial depth of shock pressure decay zone
!
!   Output
!   * xsfints   -- X position of intersection point between shock zone
!   and stream tube
! 
!***

!**********************************************************************************************************
!
!  Unit Name   : regolith_shock_damage_zone(user,crater,rm,eradi,depthb,xsfints)
!  Unit Type   : subroutine 
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculate the intersection point between shock damage
!  zone and the stream
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : xsfints  ::  
!           
! 
!  Notes       :  
!
!***********************************************************************************************************
subroutine regolith_shock_damage_zone(crater,rm,eradi,depthb,xsfints)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_shock_damage_zone
   implicit none

   ! Arguments
   type(cratertype),intent(in) :: crater
   real(DP),intent(in)         :: rm, eradi, depthb
   real(DP),intent(out)        :: xsfints
   
   ! Parameters for shock pressure calculation based on numerical modeling
   ! results (Montuex and Pierazo)
   real(DP),parameter          :: rho_dunite = 3320.0
   real(DP),parameter          :: c_dunite   = 6500.0
   real(DP),parameter          :: s_dunite   = 0.9
   real(DP),parameter          :: n          = -2.85
   real(DP)                    :: pmax, up
   real(DP)                    :: rsh, xshints0
   real(DP)                    :: q1, q2, q3, thetaq 

   ! Calculate the radius of shock pressure decay zone by a given shock 
   ! P(r) = Pmax * (r/r_p)**(n), where n is negative and we use -3 (fast decay
   ! that may be accounted for by acoustic fluidization).
   ! Make Pf equal to P(r) and solve the variable of "r", which becomes "rsh" in
   ! the following. We assume the projectile and target are the same materials.
   ! The maximum peak shock pressure is estimated from planar impact 
   up       = crater%impvel * 0.5 
   pmax     = rho_dunite * (c_dunite + s_dunite * up) * up
   rsh      = (crater%imp/2.0) * (PF / pmax)**(1.0/n)
   xshints0 = sqrt(rsh**2 - (crater%imp / 2.0)**2)

   if (eradi <= xshints0) then
      xsfints = xshints0
   else if (eradi > xshints0) then
           q1      =  1.0 / (1.0 + 2.0 * depthb / eradi)
           q2      = -1.0 - q1
           q3      = (1.0 + (depthb**2 - rsh**2)/eradi**2) * q1
           thetaq  = acos( -0.5 * q2 - 0.5 * sqrt(q2**2 - 4.0 * q3) )
           xsfints = eradi * (1.0 - cos(thetaq)) * sin(thetaq)
   end if

   return

end subroutine regolith_shock_damage_zone
