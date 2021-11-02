!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_blanket
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Scales crater ejecta blanket
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
!           
! 
!  Notes       : Ballistic trajectory equations from Fundamentals of Astrodynamics by Bate, Mueller, and White (1971). Sec. 6.2  
!
!**********************************************************************************************************************************
subroutine ejecta_blanket(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
   use module_globals
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_blanket
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: erad 
   real(DP),intent(out) :: lrad,vejsq,ejang
   logical,intent(inout) :: firstrun

   ! Internal variables
   real(DP),save :: ytrans,ytrans1,ytrans2,yterm
   real(DP),save :: ce1sq,ce2sq,ce2asq,tstr,trho,mu
   real(DP),save :: PItwo,PIfour,pvelv,p
   real(DP) :: gterm
   real(DP) :: vej1sq
   !real(DP) :: Q,cos2phi,Psi


   ! Executable code


   ! Compute some auxiliary quantities
   if (firstrun) then
      select case(crater%strflag)
      case(0)
         tstr = user%ybar_r
         trho = user%trho_r
         mu = user%mu_r
      case(1)
         tstr = user%ybar_b
         trho = user%trho_b
         mu = user%mu_b
      end select


      ! find velocity equation constant terms
      pvelv = crater%impvel * crater%sinimpang
      ce1sq = (2._DP/(ct*ct)) * mu**2 / (1._DP + 2 * mu + mu**2)
      PItwo = (user%gaccel * crater%imp * 0.5_DP) / (pvelv**2)
      PIfour = trho / user%prho
      ytrans1 = trho * pvelv**2
      ytrans2 = PItwo * (PIfour**(-THIRD))
      ytrans = ytrans1 * (ytrans2**(2._DP / (mu + 2._DP)))
      ce2asq = ce1sq * abs((user%gaccel * crater%grad * trho) / (tstr + ytrans + vsmall))
      p = 6.0_DP / mu
      ce2sq = ce2asq * (crater%grad / crater%rad)**p
      yterm = ce2sq * (tstr / trho)
      firstrun = .false.

   end if

   ! define ejection velocity at radius (erad)
   vej1sq = max(ce1sq * abs(user%gaccel * crater%grad) * (crater%grad / erad)**p,0._DP)
   gterm = ce1sq * (user%gaccel * erad)
   vejsq = max(vej1sq - gterm - yterm,0._DP)

   ! Ejection angle
   ejang=DEG2RAD * (55._DP - 20 * (erad / crater%grad))


   ! compute spherical surface ejecta landing distance for (erad)
   ! Q = vejsq / (user%gaccel * user%trad)              
   ! cos2phi = cos(ejang)**2
   ! Psi = (1._DP-Q*cos2phi)/sqrt(1._DP+Q*(Q-2._DP)*cos2phi)
   ! Psi = 2*acos(Psi)             ! Free-flight range angle
   ! lrad = erad + Psi * user%trad

   ! flat plane landing distance calculation
   lrad = erad + vejsq * sin(2 * ejang) / user%gaccel

   return
end subroutine ejecta_blanket

