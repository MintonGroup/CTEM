!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_distance_estimate
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Estimates the ejecta distance quickly to determine whether it is worth computing the whole table
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
subroutine ejecta_distance_estimate(user,crater,domain,ejdis_estimate)
   use module_globals
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_distance_estimate
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in)    :: domain
   real(DP),intent(out) :: ejdis_estimate

   ! Internal variables
   real(DP),parameter :: ERAD_DELTA_FAC = 0.05_DP
   real(DP) :: erad1,erad2,erad3,lrad1,lrad2,lrad3,thick1,thick2,erad_delta,vejsq,ejang,slope
   logical :: firstrun = .true.

   ! Executable code
   erad_delta = ERAD_DELTA_FAC * crater%rad
   erad1 = crater%rad
   vejsq = 0.0_DP
   thick1 = -1.0
   call ejecta_blanket(user,crater,domain,erad1,lrad1,vejsq,ejang,firstrun)
   do while (erad1 > 0.0_DP)
      erad2 = erad1 - erad_delta
      call ejecta_blanket(user,crater,domain,erad2,lrad2,vejsq,ejang,firstrun)
      call ejecta_thickness(user,crater,erad1,erad2,lrad1,lrad2,thick2)
      if (thick2 < domain%small) exit
      erad1 = erad2 
      lrad1 = lrad2
      thick1 = thick2
   end do
   if (thick1 < 0._DP) then
      thick1 = thick2
      erad3 = erad2 - erad_delta

      ! Get distance at three points
      call ejecta_blanket(user,crater,domain,erad2,lrad2,vejsq,ejang,firstrun)

      ! Get thickness at two locations
      call ejecta_thickness(user,crater,erad2,erad3,lrad2,lrad3,thick2)
   end if

   ! Extrapolate in logspace out to the minimum thickness
   slope = log(lrad2 / lrad1) / log(thick2 / thick1)
   ejdis_estimate = exp(slope * log(domain%small / thick1) + log(lrad1))
   ejdis_estimate = min(ejdis_estimate, crater%fcrat * user%ejecta_truncation)

   return
end subroutine ejecta_distance_estimate
