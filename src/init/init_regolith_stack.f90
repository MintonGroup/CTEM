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
subroutine init_regolith_stack(user,surf,popflag)
   use module_globals
   use module_regolith
   use module_init, EXCEPT_THIS_ONE => init_regolith_stack
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   INTEGER(I4B),DIMENSION(:,:),INTENT(INOUT)   :: popflag
   type(regodatatype) :: mare,highland
   integer(I4B) :: k,xp,yp,maresize

   ! Internal variables
   INTEGER(I4B) :: allocstat

   !call init_regolith_parab(user,surf)
   !=======================================
   ! Initialize the grid space  
   !=======================================
   DO yp = 1, user%gridsize
      DO xp = 1, user%gridsize

         IF (.NOT. ASSOCIATED(surf(xp,yp)%regolayer)) THEN
            ALLOCATE(surf(xp,yp)%regolayer, STAT=allocstat)
            IF (allocstat == 0) THEN
               NULLIFY(surf(xp,yp)%regolayer%next)

               IF (xp <= user%gridsize/2) THEN
                  highland%thickness = 1000.0_DP
                  highland%meltfrac  = 0._DP 
                  highland%comp      = 0._DP
               ELSE
                  highland%thickness = 5000.0_DP
                  highland%meltfrac  = 0._DP
                  highland%comp      = 0.0_DP
               END IF 

               surf(xp,yp)%regolayer%regodata = highland
            ELSE
               WRITE(*,*) 'Exhausted memory.'
            END IF
         ELSE
            WRITE(*,*) 'Initialization went wrong ...'
         END IF

      END DO
   END DO

   DO yp = 1,user%gridsize
      DO xp = 1,user%gridsize
         IF (xp <= user%gridsize/2) THEN
            mare%thickness = 4000.0_DP
            mare%meltfrac  = 0._DP
            mare%comp      = 1.0_DP
            call regolith_push(surf(xp,yp),mare,popflag(xp,yp))
         END IF
      END DO
   END DO

   !do yp = 1,user%gridsize
   !   do xp = 1,user%gridsize
   !         highland%thickness = 1000.0_DP
   !         highland%meltfrac  = 0._DP 
   !         highland%comp      = 0._DP
   !         call regolith_push(surf(xp,yp),highland)
   !         do k = 1, 5
   !         highland%thickness = 10.0_DP
   !         highland%meltfrac  = 0._DP
   !         highland%comp      = 0.0_DP
   !         call regolith_push(surf(xp,yp),highland)
   !         end do
   !         do k = 1, 10
   !         mare%thickness = 1.0_DP
   !         mare%meltfrac  = 0._DP
   !         mare%comp      = 1.0_DP
   !         call regolith_push(surf(xp,yp),mare)
   !         end do
   !   end do
   !end do

   return
end subroutine init_regolith_stack

