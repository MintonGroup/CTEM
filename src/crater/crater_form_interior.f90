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
subroutine crater_form_interior(user,surfi,crater,lradsq,newelev,melev)
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
   real(DP),intent(in) :: newelev,melev

   ! Internal variables
   real(DP) :: cform,newdem,elchange,pikeD
   integer(I4B) :: layer

   ! A list for poped data 
   type(regolisttype),pointer :: poppedlist

   ! Executable code

   !change digital elevation map
   cform = crater%vcorr - (crater%parab * lradsq)
   newdem = newelev - cform

   pikeD = 1.044e3_DP * (crater%fcrat * 1e-3_DP)**(0.301_DP) ! Pike (1977)
   !write(*,*) melev,pikeD
   if ((crater%fcrat > crater%cxtran * 2) .and. newdem < (melev - pikeD)) then
      newdem = melev - pikeD ! Flatten out the bottom of the crater
   end if
   if (newdem < (melev - user%deplimit)) then
      newdem = melev - user%deplimit ! Flatten out the bottom of the crater
      do layer = 1,user%numlayers ! Remove all pre-existing craters from this current pixel
         call util_remove_from_layer(surfi,layer)
      end do
   end if
   elchange  = newdem - surfi%dem
   surfi%dem = newdem
   !write(*,*) newdem
   !read(*,*)
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

