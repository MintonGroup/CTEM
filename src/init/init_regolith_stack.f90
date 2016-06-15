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
   integer(I4B) :: k,xp,yp,maresize

   ! Internal variables
   integer(I4B) :: allocstat

   !call init_regolith_parab(user,surf)
   !=======================================
   ! Initialize the grid space  
   !=======================================
   do yp = 1, user%gridsize
      do xp = 1, user%gridsize

         if (.not. associated(surf(xp,yp)%regolayer)) then
            allocate(surf(xp,yp)%regolayer, STAT=allocstat)
            if (allocstat == 0) then
               nullify(surf(xp,yp)%regolayer%next)
               bedrock%thickness = VBIG
               bedrock%meltfrac  = 0._DP 
               bedrock%comp      = 0._DP
               surf(xp,yp)%regolayer%regodata = bedrock
            else
               write(*,*) 'Exhausted memory.'
            end if
         else
            write(*,*) 'Initialization went wrong ...'
         end if

      end do
   end do

   return
end subroutine init_regolith_stack

