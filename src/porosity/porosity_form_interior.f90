!**********************************************************************************************************************************
!
!  Unit Name   : porosity_form_interior
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Compute the thickness of porosity development inside the crater. 
!       
!  Input
!    Arguments : user, surfi, crater, elchange, lradsq, lradsq 
!              : lradsq   = the distance squared from the crater center (defined in crater_emplace.f90).
!              
!  Output
!    Arguments : surfi
!           
!  Notes       : The thickness of ejecta coverage is added to surfi%ejcov 
!              : The current version only considers a binary mode: one being porosity, the other being non-porosity
!
!  Developer   : Toshi Hirabayashi (12/08/16)
!**********************************************************************************************************************************
subroutine porosity_form_interior(user, surfi, crater, lradsq)
   use module_globals
   use module_util
   use module_porosity, EXCEPT_THIS_ONE => porosity_form_interior
   implicit none

   ! Arguments
   type(usertype),intent(in)    :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in)  :: crater
   real(DP),intent(in)          :: lradsq

   ! Internal variables
   type(regodatatype)          :: nlayer  ! Dummy parameters to search for the layers.  
   type(regolisttype), pointer :: current => null() ! Dummy parameters to search for the layers. 

   ! Internal variables
   real(DP) :: newdem ! Same as surfi%dem
   ! ht: the depth of a transient crater
   ! parabht: a parabolic factor for ht.   
   ! depthht: the depth of a transient crater at some point
   real(DP) :: ht, parabht, depthht, newtrn
 
   ! Executable code
   ! Interface
   newdem = surfi%dem

   ! compute Transient crater depth
   ht = TRNRATIO * crater%fcrat      !The depth of a transient crater   
   parabht = ht / ((crater%frad)**2) !The factor of a parabolic profile  
   depthht = ht - (parabht * lradsq) !The depth at a given pixel.  
   newtrn  = crater%melev - depthht  !The elevation of the bottom of the transient crater at a given pixel.

	nlayer%porosity = 0.2;            !The porosity of the first value. This value is an assumed value. 
	nlayer%depth    = newtrn; 

	! Currently consider only the first layer. 
	if (surfi%porolayer%regodata%porosity == nlayer%porosity) then
		! just update if the newly calculated layer is deeper than the stored one. 
		if (surfi%porolayer%regodata%depth > nlayer%depth) then
			surfi%porolayer%regodata%depth    = nlayer%depth
			surfi%porolayer%regodata%porosity = nlayer%porosity
		end if
	else
		! Linked list if there is only one porosity layer
		call util_push(surfi%porolayer, nlayer)
	end if 
	
   return
end subroutine porosity_form_interior

