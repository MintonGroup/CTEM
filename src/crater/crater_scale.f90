!**********************************************************************************************************************************
!
!  Unit Name   : crater_scale
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Scales crater from projectile size
!
!  Input
!    Arguments : psize   : Projectile size (m)
!              : strflag : Strength flag (0 = regolith; 1 = bedrock)
!              : sinimpang  : Sine of the impact angle
!              : impvel  : Impact velocity (m/s)
!
!  Output
!    Arguments : crad  : Crater radius (m)
!              : grad  : Gravity scaled radius (m)
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine crater_scale(user,psize,crad,grad,strflag,sinimpang,impvel)
   use module_globals
   use module_crater, EXCEPT_THIS_ONE => crater_scale
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   real(DP),intent(in) :: psize,sinimpang,impvel
   real(DP),intent(out) :: crad,grad
   integer(I4B),intent(in) :: strflag

   ! Internal variables
   real(DP)           :: cvol,cvolg
   real(DP)           :: pivol,pivolg,pitwo,pithree,pifour
   real(DP)           :: prad,pvelv,pmass,tstr,trho,mu,kv
   real(DP)           :: c1,c2
   !real(DP)           ::  ct,exdep,exfac,massej,tform

   ! Executable code

   ! set target strength
   select case(strflag)
   case(0)
      tstr = user%ybar_r
      trho = user%trho_r
      mu = user%mu_r
      kv = user%kv_r
   case(1)
      tstr = user%ybar_b
      trho = user%trho_b
      mu = user%mu_b
      kv = user%kv_b
   case default
      tstr = user%ybar_r
      trho = user%trho_r
      mu = user%mu_r
      kv = user%kv_r
   end select

   ! Compute some auxiliary quantites
   pvelv = impvel * sinimpang
   prad = 0.5_DP * psize
   pmass = 4 * THIRD * pi * user%prho * prad**3
   c1 = 1._DP + 0.5_DP * mu
   c2 = (-3 * mu)/(2._DP + mu)

   ! Find dimensionless quantities
   pitwo = (user%gaccel * prad)/(pvelv**2)
   pithree = tstr / (trho * (pvelv**2))
   pifour = trho / user%prho
   pivol = kv * ((pitwo * (pifour**(-THIRD))) + (pithree**c1))**c2
   pivolg = kv * (pitwo * (pifour**(-THIRD)))**c2
   
   ! find transient crater volume and radii (depth = 1/3 diameter)
   cvol = pivol * (pmass / trho)
   cvolg = pivolg * (pmass / trho)
   crad = (3 * cvol / PI)**(THIRD)
   grad = (3 * cvolg / PI)**(THIRD)

   ! find the crater formation time
   !tform = CT * ((crad / grad)**((mu + 1._DP) / mu)) * sqrt(abs(grad / user%gaccel))

   ! find the excavation parameters
   !exdep = (2 * crad) / exfac
   !massej = trho * 0.5_DP * PI * exdep * (crad**2)

   return
end subroutine crater_scale

