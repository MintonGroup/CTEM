!**********************************************************************************************************************************
!
!  Unit Name   : crater_form_interior
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
subroutine crater_form_interior(user,surfi,crater,x_relative, y_relative ,newelev,deltaMi)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_form_interior
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(inout) :: crater
   real(DP),intent(in) :: x_relative, y_relative 
   real(DP),intent(in) :: newelev
   real(DP),intent(out) :: deltaMi

   ! Internal variables
   real(DP) :: cform,newdem,elchange,r
   integer(I4B) :: layer

   ! A list for popped data 
   type(regolisttype),pointer :: poppedlist


   ! Empirical crater shape parameters from Fassett et al. (2014)
   real(DP),parameter :: r_floor = 0.2_DP
   real(DP),parameter :: simple_depth_diam = 0.181_DP
   real(DP),parameter :: r_rim = 0.98_DP
   real(DP),parameter :: inner_c0 = -0.229_DP 
   real(DP),parameter :: inner_c1 =  0.228_DP 
   real(DP),parameter :: inner_c2 =  0.083_DP 
   real(DP),parameter :: inner_c3 = -0.039_DP

   real(DP),parameter :: outer_c0 =  0.188_DP
   real(DP),parameter :: outer_c1 = -0.187_DP
   real(DP),parameter :: outer_c2 =  0.018_DP 
   real(DP),parameter :: outer_c3 =  0.015_DP

   real(DP) :: c0,c1,c2,c3,flrad,rh,fld


   ! Executable code
   r = sqrt(x_relative**2+y_relative**2) / crater%frad

   rh = crater%rimheight 
   fld = -crater%floordepth 
   flrad = 0.5_DP * crater%floordiam / crater%frad 

   
   ! Use polynomial crater profile similar to that of Fassett et al. (2014), but the parameters are set by the crater dimensions
   c1 = (fld - rh) / (flrad + flrad**2 / 3._DP - flrad**3 / 6._DP - 7._DP / 6._DP)
   c0 = rh - (7._DP / 6._DP) * c1
   c2 = c1 / 3._DP
   c3 = -c2 / 2._DP

   if (r < flrad) then
      cform = fld 
   else
      cform = c0 + c1 * r + c2 * r**2 + c3 * r**3 
   end if
   ! TEMP UNTIL I WRITE THE SOLUTION PROPERLY
   if (cform > 0.0_DP .and. r * crater%frad < crater%ejrad) crater%ejrad = r * crater%frad

   !if (r < r_floor) then
   !   cform = -simple_depth_diam  * crater%fcrat
   !else if (r < r_rim) then
   !   cform =  (inner_c0 + inner_c1 * r + inner_c2 * r**2 + inner_c3 * r**3) * crater%fcrat 
   !else 
   !   cform =  (outer_c0 + outer_c1 * r + outer_c2 * r**2 + outer_c3 * r**3) * crater%fcrat
   !end if      
   newdem = newelev + cform 

   !if (crater%fcrat > crater%cxtran * 2) then 
   if (crater%morphtype == "COMPLEX") then
      ! Make this crater complex

      !if (r < r_rim) then
         if (newdem < (crater%melev - crater%floordepth)) then ! Flatten out the bottom of the crater
            do layer = 1,user%numlayers ! Remove all pre-existing craters from this current pixel
               call util_remove_from_layer(surfi,layer)
            end do
            newdem = crater%melev - crater%floordepth

            
         end if
      !end if

   end if


   newdem = min(newdem,surfi%dem) ! Only allow excavation, no deposition
   elchange  = newdem - surfi%dem
   deltaMi = elchange
   surfi%dem = newdem
   
   !change ejecta coverage
   surfi%ejcov = max(surfi%ejcov + elchange,0.0_DP)

   if (user%doregotrack) then
      call util_traverse_pop(surfi%regolayer,abs(elchange),poppedlist)
      call util_destroy_list(poppedlist)
   end if


   return
end subroutine crater_form_interior

