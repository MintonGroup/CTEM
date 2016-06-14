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
   type(regodatatype) :: mare,highland
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

               if (xp <= user%gridsize/2) then
                  highland%thickness = 1000.0_DP
                  highland%meltfrac  = 0._DP 
                  highland%comp      = 0._DP
               else
                  highland%thickness = 5000.0_DP
                  highland%meltfrac  = 0._DP
                  highland%comp      = 0.0_DP
               end if 

               surf(xp,yp)%regolayer%regodata = highland
            else
               write(*,*) 'Exhausted memory.'
            end if
         else
            write(*,*) 'Initialization went wrong ...'
         end if

      end do
   end do

   do yp = 1,user%gridsize
      do xp = 1,user%gridsize
         if (xp <= user%gridsize/2) then
            mare%thickness = 4000.0_DP
            mare%meltfrac  = 0._DP
            mare%comp      = 1.0_DP
            call util_push(surf(xp,yp)%regolayer,mare)
         end if
      end do
   end do

   return
end subroutine init_regolith_stack

