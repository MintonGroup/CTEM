!**********************************************************************************************************************************
!
!  Unit Name   : crater_find_visible
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds the visible crater parabolic parameters, rim, and rim upturn distance
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
subroutine crater_find_visible(user,crater,domain)
   use module_globals
   use module_crater, EXCEPT_THIS_ONE => crater_find_visible
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in)    :: domain

   ! Internal variables
   real(DP) :: cdepth,lrad,cform

   ! Executable code
   cdepth = DDRATIO * crater%fcrat
   if (crater%fcrat <= crater%cxtran) then
      crater%rheight = RDRATIO * crater%fcrat
   else
      crater%rheight = (RDRATIO*crater%cxtran)+(RDRATIO*((crater%fcrat-crater%cxtran)**(0.35_DP)))
   endif

   crater%vcorr  = cdepth - crater%rheight
   crater%parab  = cdepth / ((crater%frad)**2)

   !find rim for counting purposes
   crater%frim = RIMFAC * crater%frad

   ! find rim upturn distance
   lrad = crater%frad
   crater%rimdis = domain%side
   do while (lrad <= domain%side)
      cform = crater%rheight * ((crater%frad/lrad)**RIMDROP)
      if (cform < domain%small) then
         crater%rimdis = lrad
         exit
      end if
      lrad = lrad + (user%pix * SUBPIXFAC)
   end do
   ! Get pixel space values
   crater%frimpx = int(crater%frim/user%pix) + 1
   crater%rimdispx = int(crater%rimdis/user%pix)  + 1

   return
end subroutine crater_find_visible

