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
function regolith_shock_damage(erad,deltar,xmints,xsfints,xleft,xright) result(vsh)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_shock_damage
   implicit none

   ! Arguments
   real(DP),intent(in)         :: erad, deltar, xmints, xsfints, xleft, xright
   real(DP)                    :: vsh
   real(DP)                    :: x_low_sh, x_up_sh  
 
   if (xsfints > xmints) then
      x_low_sh         = max(xmints,xleft)
      x_up_sh          = min(xsfints,xright)
      vsh              = regolith_streamtube_volume_func(erad,x_low_sh,x_up_sh,deltar)
   else
      vsh              = 0.0_DP
   end if

   return

end function regolith_shock_damage
