!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_thickness
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Scales crater ejecta thickness
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
subroutine ejecta_thickness(user,crater,erad1,erad2,lrad1,lrad2,thick)
   use module_globals
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_thickness
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: erad1,erad2,lrad1,lrad2
   real(DP),intent(out) :: thick

   ! Internal variables
   real(DP) :: excav,landarea,ce1sq,edDt,mu
   !real(DP) :: Psi1,Psi2,h1,h2,a1,a2

   ! Executable code
   select case(crater%strflag)
   case(0)
      mu = user%mu_r
   case(1)
      mu = user%mu_b
   case default
      mu = user%mu_r
   end select

   if (crater%imp < user%basinimp) then
      !  compute ejecta thickness 
      ce1sq = (2._DP / (CT * CT)) * mu**2 / (1._DP + 2 * mu + mu**2)
      excav = abs(0.5_DP * ce1sq * PI * (erad2**3 - erad1**3))
   else
      edDt = 1.150_DP ! Calibrated to Orientale using Fasset et al. 2011
      excav = abs(0.5_DP * edDt * PI * (erad2**3 - erad1**3))
   end if


   ! Compute area as differences between areas of two spherical caps
   !Psi1=(lrad1-erad1)/user%trad
   !Psi2=(lrad2-erad2)/user%trad

   !a1=user%trad*sin(Psi1)
   !a2=user%trad*sin(Psi2)
   !if (Psi1 > 0.5_DP*PI) then
   !   h1 = user%trad-sqrt(user%trad**2-a1**2)
   !else
   !   h1 = user%trad+sqrt(user%trad**2-a1**2)
   !end if
   !if (Psi2 > 0.5_DP*PI) then 
   !   h2 = user%trad-sqrt(user%trad**2-a2**2)
   !else
   !   h2 = user%trad+sqrt(user%trad**2-a2**2)
   !end if
   !landarea = abs(2*PI*user%trad*(h1-h2))

   landarea = PI * (lrad2*lrad2 - lrad1*lrad1)

   thick = excav / (landarea+VSMALL)
   if (thick < VSMALL) thick = VSMALL

   return
end subroutine ejecta_thickness

