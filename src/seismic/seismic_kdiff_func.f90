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
function seismic_kdiff_func(user,crater,var,gratio,firstrun,invflag) result(ans)
   use module_globals
   use module_seismic, EXCEPT_THIS_ONE => seismic_kdiff_func
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: var ! kdiff when invflag is false and lrad if invflag is true
   real(DP),intent(out) :: gratio
   logical,intent(inout) :: firstrun
   logical,intent(in) :: invflag
   real(DP) :: ans ! lrad when invflag is false and kdiff if invflag is true

   
   ! Internal variables
   real(DP) :: attenl,saccel,lrad,kdiff
   real(DP),save :: seisk,kdiffterm,saccelterm,cohaccel
   !real(DP) :: rat,lcirc,lang

   ! Seismic parameters
   real(DP),parameter :: SEISFREQ = 20.0_DP    ! seismic wave frequency
   real(DP),parameter :: SHEFF = 0.15_DP  ! seismic diffusion distance factor
   real(DP),parameter :: QFAC = 1.80_DP  ! seismic diffusion distance factor
   real(DP),parameter :: NFAC = 0.25_DP  ! seismic diffusion distance factor
   real(DP),parameter :: IMPFAC = 0.75_DP  ! seismic diffusion impactor diameter factor
   real(DP),parameter :: KSFAC = 0.30_DP  ! seismic diffusion impactor diameter factor

   ! Executable Code

   ! find acceleration
   if (firstrun) then
      seisk = THIRD * user%tvel * user%tfrac
      cohaccel = user%regcoh / user%trho_r
      kdiffterm = SHEFF * (user%seisq**QFAC) * (user%neff**NFAC) * (crater%imp**IMPFAC) * sqrt(crater%impvel) * sqrt(user%gaccel)
      kdiffterm = kdiffterm / (seisk**KSFAC)
      saccelterm = user%neff * user%prho * (crater%impvel**2) * (SEISFREQ**2) * (crater%imp**3)
      firstrun = .false.
   end if

   ! Calculate gratio
   if (invflag) then ! calculate the inverse function, without the gratio falloff factor
      kdiff = var

      !lcirc = (kdiffterm/kdiff)**(1._DP/distfac)
      !rat = lcirc/(2*user%trad)
      !if (abs(rat)<=1._DP) then
      !   lang = 2*asin(rat)
      !else
      !   lang = sign(PI,rat)
      !end if
      !lrad = user%trad * lang

      lrad = (kdiffterm / (kdiff + VSMALL))**2
      ans = lrad

      attenl = exp((-2 * SEISFREQ * (lrad**2)) / (seisk * PI * user%seisq))
      saccel = saccelterm * attenl
      saccel = (PI / SQRT3) * sqrt(abs(saccel / (user%trho_b * (lrad**2))))
      gratio = 0.1_DP * (max(saccel-cohaccel, 0._DP) / user%gaccel)
   else
      lrad = var
      kdiff = kdiffterm / (sqrt(lrad) + VSMALL)

      !lang = lrad / user%trad
      !lcirc = 2*user%trad*dsin(0.5_DP*lang)

      ! find downslope diffusion
      !kdiff = kdiffterm / (lcirc**distfac + VSMALL)

      ! Calculate gratio
      !attenl = exp((-2*SEISFREQ*(lcirc**2))/(seisk*PI*user%seisq))
      !saccel = (PI/SQRT3)*sqrt(abs(saccel/(user%trho_b*(lcirc**2))))

      attenl = exp((-2 * SEISFREQ * (lrad**2)) / (seisk * PI * user%seisq))
      saccel = saccelterm * attenl
      saccel = (PI/SQRT3) * sqrt(abs(saccel / (user%trho_b * (lrad**2))))
      gratio = 0.1_DP * (max(saccel-cohaccel, 0._DP) / user%gaccel)
      ans = kdiff*min(sqrt(gratio),1._DP)
   end if

   return
   end function seismic_kdiff_func
