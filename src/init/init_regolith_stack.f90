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
subroutine init_regolith_stack(user,surf,domain)
   use module_globals
   use module_util
   use module_init, EXCEPT_THIS_ONE => init_regolith_stack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(domaintype),intent(in)    :: domain
   integer(I4B) :: xp,yp

   ! Internal variables
   logical :: initstat

   !call init_regolith_parab(user,surf)
   !=======================================
   ! Initialize the grid space  
   !=======================================
   if (user%dothermal) then
      do yp=1,user%gridsize
         do xp=1,user%gridsize
            call util_init_array_split(user,surf(xp,yp)%regolayer,domain,initstat)
         end do
      end do
   else
      do yp=1,user%gridsize
         do xp=1,user%gridsize
            call util_init_array(user,surf(xp,yp)%regolayer,domain,initstat)
         end do
      end do
   end if

   return
end subroutine init_regolith_stack

