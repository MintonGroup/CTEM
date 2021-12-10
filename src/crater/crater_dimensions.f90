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

!   real(DP),parameter :: DDRATIO = 0.19_DP        ! Depth-diameter ratio
!   real(DP),parameter :: RDRATIO = 0.030_DP       ! Rim height to diameter ratio
   real(DP),parameter :: RIMFAC = 1.5_DP          ! Ratio of radius used for counting craters to the final rim radius

   !TODO: Add in transitional and multi-ringed basin crater morphology
   !TODO: Implement Monte Carlo model for standard errors for all dimension parameters

   ! Set the crater morphology
   select case(crater%morphtype)
      case("SIMPLE","TRANSITION") ! A hybrid model between Pike (1977) and Fassett & Thomson (2014) 
         !crater%rimheight  = 0.036_DP * (crater%fcrat * 1e-3_DP)**(1.014_DP) * 1e3_DP !* crater%fcrat ! Pike model
         crater%rimheight  = 0.043_DP * (crater%fcrat * 1e-3_DP)**(1.014_DP) * 1e3_DP !* crater%fcrat  ! Closer to Fassett & Thomson
         crater%rimwidth   = 0.257_DP * (crater%fcrat * 1e-3_DP)**(1.011_DP) * 1e3_DP !* crater%fcrat  ! Pike model
         !crater%floordepth = 0.196_DP * (crater%fcrat * 1e-3_DP)**(1.010_DP) * 1e3_DP !* crater%fcrat ! Pike model
         crater%floordepth = 0.224_DP * (crater%fcrat * 1e-3_DP)**(1.010_DP) * 1e3_DP !* crater%fcrat  ! Closer to Fassett & Thomson
         !crater%floordiam  = 0.031_DP * (crater%fcrat * 1e-3_DP)**(1.765_DP) * 1e3_DP !* crater%fcrat ! Pike model
         crater%floordiam  = 0.200_DP * (crater%fcrat * 1e-3_DP)**(1.143_DP) * 1e3_DP !* crater%fcrat  ! Fassett & Thomson for D~1km, Pike for D~20km
      case("COMPLEX","PEAKRING","MULTIRING") ! Following Pike (1977)
         crater%rimheight  = 0.236_DP * (crater%fcrat * 1e-3_DP)**(0.399_DP) * 1e3_DP !* crater%fcrat
         crater%rimwidth   = 0.467_DP * (crater%fcrat * 1e-3_DP)**(0.836_DP) * 1e3_DP !* crater%fcrat
         crater%floordepth = 1.044_DP * (crater%fcrat * 1e-3_DP)**(0.301_DP) * 1e3_DP !* crater%fcrat
         crater%floordiam  = min(0.187_DP * (crater%fcrat * 1e-3_DP)**(1.249_DP) * 1e3_DP, 0.9_DP * crater%fcrat) !* crater%fcrat
         crater%peakheight = 0.032_DP * (crater%fcrat * 1e-3_DP)**(0.900_DP) * 1e3_DP !* crater%fcrat 
   end select


   ! Executable code
   !crater%floordepth = DDRATIO * crater%fcrat
   !if (crater%fcrat <= crater%cxtran) then
   !else
   !   crater%rimheight = (RDRATIO*crater%cxtran)+(RDRATIO*((crater%fcrat-crater%cxtran)**(0.399_DP)))
   !endif

   ! Adjust the floor depth to measure from the pre-existing level surface, rather than the rim
   crater%floordepth  = crater%floordepth - crater%rimheight

   ! Calculate the radius where the inner wall meets the original pre-existing surface
   ! This is used to demark the location where excavation transitions to deposition
   crater%ejrad = crater_profile_find_r_inner_wall(user,crater) * crater%frad

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

