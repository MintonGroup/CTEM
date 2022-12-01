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
   use module_util
   use module_init, EXCEPT_THIS_ONE => init_regolith_stack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(regodatatype) :: bedrock 
   integer(I4B) :: xp,yp,k

   ! Internal variables
   logical :: initstat

   ! Temporary variable setup for initialize a pre-exising structure
   type(regodatatype) :: test_stratig

   !call init_regolith_parab(user,surf)
   !=======================================
   ! Initialize the grid space  
   !=======================================
   bedrock%thickness = user%trad
   bedrock%meltfrac  = 0._DP 
   bedrock%comp      = 0._DP
   bedrock%age(:)    = 0.0_SP

   do yp = 1, user%gridsize
      do xp = 1, user%gridsize

         !call util_init_list(surf(xp,yp)%regolayer,initstat)
         call util_init_array(surf(xp,yp)%regolayer,initstat)

         if (initstat) then
             call util_push_array(surf(xp,yp)%regolayer,bedrock)
         else
            write(*,*) 'init_regolith_stack: Initialization of regolayer failed.'
         end if

      end do
   end do

   return
end subroutine init_regolith_stack

