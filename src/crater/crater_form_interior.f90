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
subroutine crater_form_interior(user,surfi,crater,lradsq,newelev,deltaMi)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_form_interior
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: lradsq
   real(DP),intent(in) :: newelev
   real(DP),intent(out) :: deltaMi

   ! Internal variables
   real(DP) :: cform,newdem,elchange,pikeD,r,h,circrad,polymodel,f,HH
   integer(I4B) :: layer

   ! A list for poped data 
   type(regolisttype),pointer :: poppedlist

   ! Executable code

   !change digital elevation map
   r = sqrt(lradsq) / crater%frad
   ! Use empirical crater form from Fassett et al. 2014
   if (r < 0.2_DP) then
      cform = -0.181_DP * crater%fcrat
   else if (r < 0.98_DP) then
      cform =  (-0.229_DP + 0.228_DP * r + 0.083_DP * r**2 - 0.039_DP * r**3) * crater%fcrat 
   else 
      cform =  (0.188_DP - 0.187_DP * r + 0.018_DP * r**2 + 0.015_DP * r**3) * crater%fcrat
   end if      
   newdem = newelev + cform 

   pikeD = 1.044e3_DP * (crater%fcrat * 1e-3_DP)**(0.301_DP) ! Pike (1977)
   if ((crater%fcrat > crater%cxtran * 2) .and. newdem < (crater%melev - pikeD)) then
      newdem = crater%melev - pikeD ! Flatten out the bottom of the crater
   end if
   if (newdem < (crater%melev - user%deplimit)) then
      newdem = crater%melev - user%deplimit ! Flatten out the bottom of the crater
      do layer = 1,user%numlayers ! Remove all pre-existing craters from this current pixel
         call util_remove_from_layer(surfi,layer)
      end do
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

