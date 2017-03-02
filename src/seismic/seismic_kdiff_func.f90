!**********************************************************************************************************************************
!
!  Unit Name   : seismic_kdiff_func
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Seismic shaking diffusion function
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       :  
!
!**********************************************************************************************************************************
function seismic_kdiff_func(user,crater,var,gratio,invflag) result(ans)
   use module_globals
   use module_seismic, EXCEPT_THIS_ONE => seismic_kdiff_func
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(out) :: gratio
   logical,intent(in) :: invflag
   real(DP),intent(in) :: var ! kdiff when invflag is true and lrad if invflag is false

   ! returned value
   real(DP) :: ans ! lrad when invflag is true and kdiff if invflag is false
   
   ! Internal variables
   real(DP) :: attenl,saccel,lrad,kdiff
   !real(DP) :: rat,lcirc,lang

   ! Executable Code

   ! Calculate gratio
   if (invflag) then ! calculate the inverse function, without the gratio falloff
      kdiff = var
      lrad = (crater%kdiffterm / (kdiff + VSMALL))**(1._DP/DFAC)
      ans = lrad
      attenl = exp((-2 * SEISFREQ * (lrad**2)) / (user%seisk * PI * user%seisq))
      saccel = crater%saccelterm * attenl
      saccel = (PI / SQRT3) * sqrt(abs(saccel / (user%trho_b * (lrad**2))))
      gratio = 0.1_DP * (max(saccel-user%cohaccel, 0._DP) / user%gaccel)
   else
      lrad = var
      kdiff = crater%kdiffterm / ((lrad + VSMALL)**DFAC)
      attenl = exp((-2 * SEISFREQ * (lrad**2)) / (user%seisk * PI * user%seisq))
      saccel = crater%saccelterm * attenl
      saccel = (PI/SQRT3) * sqrt(abs(saccel / (user%trho_b * (lrad**2))))
      gratio = 0.1_DP * (max(saccel-user%cohaccel, 0._DP) / user%gaccel)
      ans = kdiff*min(sqrt(gratio),1._DP)
   end if

   return
   end function seismic_kdiff_func
