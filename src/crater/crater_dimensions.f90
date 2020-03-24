!**********************************************************************************************************************************
!
!  Unit Name   : crater_dimensions
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the physical dimensions of the crater
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
subroutine crater_dimensions(user,crater,domain)
   use module_globals
   use module_crater, EXCEPT_THIS_ONE => crater_dimensions
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in)    :: domain

   ! Internal variables
   real(DP) :: lrad,cform

   real(DP),parameter :: DDRATIO = 0.19_DP        ! Depth-diameter ratio
   real(DP),parameter :: RDRATIO = 0.030_DP       ! Rim height to diameter ratio
   real(DP),parameter :: RIMFAC = 1.5_DP          ! Ratio of radius used for counting craters to the final rim radius


   ! Executable code
   crater%floordepth = DDRATIO * crater%fcrat
   if (crater%fcrat <= crater%cxtran) then
      crater%rimheight = RDRATIO * crater%fcrat
   else
      crater%rimheight = (RDRATIO*crater%cxtran)+(RDRATIO*((crater%fcrat-crater%cxtran)**(0.399_DP)))
   endif

   crater%vcorr  = crater%floordepth - crater%rimheight
   crater%parab  = crater%floordepth / ((crater%frad)**2)

   !find rim for counting purposes
   crater%frim = RIMFAC * crater%frad

   ! find rim upturn distance
   lrad = crater%frad
   crater%rimdis = domain%side
   do while (lrad <= domain%side)
      cform = crater%rimheight * ((crater%frad/lrad)**RIMDROP)
      if (cform < domain%small) then
         crater%rimdis = lrad
         exit
      end if
      lrad = lrad + (user%pix * SUBPIXFAC)
   end do
   ! Get pixel space values
   crater%fradpx = int(crater%frad/user%pix) + 1
   crater%rimdispx = int(crater%rimdis/user%pix)  + 1

   return
end subroutine crater_dimensions

