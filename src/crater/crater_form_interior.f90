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
   use module_porosity
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
   !cform = crater%vcorr - (crater%parab * lradsq)
   !if (lradsq < crater%rad**2) then
   r = sqrt(lradsq) / crater%frad
   h = 2.5_DP*crater%rad !0.5_DP*crater%fcrat
   circrad = 0.5_DP * h + (2*crater%rad)**2 / (8 * h)
   cform = sqrt(circrad**2 - lradsq) + (circrad - h) + 0.15_DP*crater%frad
                  
   newdem = newelev - cform

   pikeD = 1.044e3_DP * (crater%fcrat * 1e-3_DP)**(0.301_DP) ! Pike (1977)
   !write(*,*) melev,pikeD
   if ((crater%fcrat > crater%cxtran * 2) .and. newdem < (crater%melev - pikeD)) then
      newdem = crater%melev - pikeD ! Flatten out the bottom of the crater
   end if
   if (newdem < (crater%melev - user%deplimit)) then
      newdem = crater%melev - user%deplimit ! Flatten out the bottom of the crater
      do layer = 1,user%numlayers ! Remove all pre-existing craters from this current pixel
         call util_remove_from_layer(surfi,layer)
      end do
   end if
   elchange  = newdem - surfi%dem
   deltaMi = elchange
   surfi%dem = newdem
   
  if (user%doporosity) then
      call porosity_form_interior(user,surfi,crater,elchange,lradsq,newelev)
   else
      !change ejecta coverage
      surfi%ejcov = max(surfi%ejcov + elchange,0.0_DP)
   end if    
     
   if (user%doregotrack) then
      call util_traverse_pop(surfi%regolayer,abs(elchange),poppedlist)
      call util_destroy_list(poppedlist)
   end if


   return
end subroutine crater_form_interior

