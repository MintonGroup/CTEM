!**********************************************************************************************************************************
!
!  Unit Name   : init_porosity_stack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the porosity stack to null
!  
!
!  Input
!    Arguments : porolayer
!
!  Output
!    Arguments : porolayer
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine init_porosity_stack(user,surf)
   use module_globals
   use module_util
   use module_init, EXCEPT_THIS_ONE => init_porosity_stack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(regodatatype) :: bedrock    
   integer(I4B) :: xp, yp

   ! Internal variables
   logical :: initstat

   do yp = 1, user%gridsize
      do xp = 1, user%gridsize

      	call util_init_list(surf(xp,yp)%porolayer, initstat)
         
      end do
   end do

   return
end subroutine init_porosity_stack

