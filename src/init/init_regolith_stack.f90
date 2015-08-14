!**********************************************************************************************************************************
!
!  Unit Name   : init_regolith_stack
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the regolith stack to null
!  
!
!  Input
!    Arguments : regolayer
!
!  Output
!    Arguments : regolayer
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine init_regolith_stack(user,surf)
   use module_globals
   use module_regolith
   use module_init, EXCEPT_THIS_ONE => init_regolith_stack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(regolayertype),pointer :: topregolayer
   type(regolayertype) :: cleanrego

   ! Internal variables
   integer(I4B) :: xp,yp

   cleanrego%thickness = huge(0._DP)
   cleanrego%mixfrac = 0._DP

!   do xp = 1,user%pix
!      do yp = 1,user%pix
!         nullify(surf(xp,yp)%regolayer)
!         allocate(surf(xp,yp)%regolayer)
!         nullify(surf(xp,yp)%regolayer%next)
!         call regolith_push(surf(xp,yp)%regolayer,cleanrego)
!      end do
!   end do

   return
end subroutine init_regolith_stack

