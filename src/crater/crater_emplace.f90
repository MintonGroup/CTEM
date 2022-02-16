!****f* crater/crater_emplace
! Name
!   crater_emplace -- 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   * module_crater
!   
!   call crater_emplace(user,surf,crater,domain,deltaMtot)
!
! DESCRIPTION
!    
! 
! ARGUMENTS
!   Input
!   * user    -- User input parameters
!   * surf    -- Surface grid
!   * crater  -- Crater dimension container
!   * domain  -- Simulation domain variable container
!   * deltaMtot -- Total displaced mass used to calculate mass-conserving ejecta blanket
!
!   Output
!   * surf   -- Outputs the new ejecta blanket onto the grid
!   * crater -- May affects the value of the maximum affected distance
!   * domain -- 
!
!   teswelktjrlkgjdlr
! 
! Notes
!
!***

!**********************************************************************************************************************************
!
!  Unit Name   : crater_emplace
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds the visible crater parabolic parameters, rim, and  rim upturn distance
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
subroutine crater_emplace(user,surf,crater,domain,deltaMtot)
   use module_globals
   use module_util
   use module_porosity   
   use module_crater, EXCEPT_THIS_ONE => crater_emplace
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(inout) :: domain
   real(DP),intent(out) :: deltaMtot

   ! Internal variables
   real(DP) :: lradsq,newelev, x_relative, y_relative 
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq,deltaMi,rimheight
   logical :: lastloop

   ! Executable code

   ! determine area to effect
   ! First make the interior of the crater
   inc = max(min(crater%fradpx,PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   deltaMtot = 0.0_DP !ejbmass
   incsq = inc**2
   ! This loop may not be parallelizable because of the linked list operation inside crater_form_interior
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         iradsq = i**2 + j**2 
         ! find distance from crater center
         
         ! find elevation and grid point
         newelev = crater%melev + ((i * crater%xslp) + (j * crater%yslp)) * user%pix
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j

         xp = xpi * user%pix
         yp = ypi * user%pix
         
         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)
         x_relative = (crater%xl - xp)
         y_relative = (crater%yl - yp)
         
         lradsq = x_relative**2 + y_relative**2

         if (lradsq > crater%frad**2) cycle
         call crater_form_interior(user,surf(xpi,ypi),crater,x_relative, y_relative,newelev,deltaMi)
         deltaMtot = deltaMtot + deltaMi

         ! do porosity computation if (user%doporosity)
         ! It is still important to consider the physical meaning of frad and rad. 
         ! frad is the final crater, while rad is the transient crater. 
         ! which one should be reasonable here. 
         if (user%doporosity) then
            call porosity_form_interior(user, surf(xpi,ypi), crater, lradsq)
         end if
     
      end do
   end do

   domain%tallycoverage = domain%tallycoverage + int(fradsq * PI / user%pix**2)
   domain%subpixelcoverage = domain%subpixelcoverage + int(fradsq * PI / user%pix**2)

   return
end subroutine crater_emplace

