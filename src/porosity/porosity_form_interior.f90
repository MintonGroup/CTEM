!**********************************************************************************************************************************
!
!  Unit Name   : porosity_form_interior
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Compute the thickness of ejecta coverage due to a transient crater at each pixel. 
!       
!
!  Input
!    Arguments : user, surfi, crater, elchange, lradsq, lradsq 
!              : elchange = the elevation change from the previous run to the latest one (defined in crater_form_interior.f90). 
!              : lradsq   = the distance squared from the crater center (defined in crater_emplace.f90).
!              : newelev  = the new elevation in the previous run at a given pixel (defined in crater_emplace.f90). 
!              
!  Output
!    Arguments : surfi
!           
! 
!  Notes       : The thickness of ejecta coverage is added to surfi%ejcov 
!              : The current version only considers a binary mode: one being porosity, the other being non-porosity
!
!  Developer   : Toshi Hirabayashi (01/16/16)
!**********************************************************************************************************************************
subroutine porosity_form_interior(user,surfi,crater,elchange,lradsq,newelev)
   use module_globals
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: lradsq,elchange,newelev

   ! Internal variables
   real(DP) :: newdem ! Same as surfi%dem
   ! ht: the depth of a transient crater
   ! parabht: a parabolic factor for ht.   
   ! depthht: the depth of a transient crater at some point
   real(DP) :: ht, parabht, depthht, newtrn
 
   ! Executable code
   ! Interface
   newdem=surfi%dem

   !compute Transition crater depth
   ht=0.3*crater%fcrat               !The depth of a transient crater   
   parabht = ht / ((crater%frad)**2) !The factor of a parabolic profile  
   depthht = ht - (parabht * lradsq) !The depth at a given pixel.  
   newtrn  = newelev - depthht       !The elevation of the bottom of the transient crater at a given pixel.

   ! Impact makes regolith
   if (surfi%ejcov < max(newdem - newtrn,0.0_DP)) then !This condition was made because depthht is the magnitude
                                            !while newdem is the elevation. Thus, to have the difference between them, 
                                            !depthht + newdem. 
   		surfi%ejcov = max(surfi%ejcov + newdem - newtrn,0.0_DP)
   else 
   		surfi%ejcov = max(surfi%ejcov + min(elchange,0.0_DP),0.0_DP) 
   end if 	

   return
end subroutine porosity_form_interior

