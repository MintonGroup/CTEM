!**********************************************************************************************************************************
!
!  Unit Name   : init_surf
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the surface arrays
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine init_surf(user,surf)
   use module_globals
   use module_init, EXCEPT_THIS_ONE => init_surf
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(out) :: surf

   ! Internal variables
   integer(I4B) :: layer

   surf%ejcov  = 0.0_DP
   surf%dem    = 0.0_DP
   do layer = 1,user%numlayers
      surf%diam(layer)   = 0.0_DP
      surf%xl(layer)     = 0.0_SP
      surf%yl(layer)     = 0.0_SP
   end do
   !if (user%docrustal_thinning) surf%mantle = 0._DP

   if (user%doregotrack) call init_regolith_stack(user,surf)
   
   ! If doporosity, call init_porosity_stack to define the porolayer linked list. 
 	if (user%doporosity)  call init_porosity_stack(user,surf)

   return
end subroutine init_surf
