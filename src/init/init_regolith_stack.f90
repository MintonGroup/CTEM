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
   type(regolayertype) :: mare,highland
   integer(I4B) :: i,xp,yp,maresize

!   cleanrego%thickness = 1500.0_DP !5000.0_DP
!   cleanrego%meltfrac = 0._DP
!   cleanrego%comp = 0._DP
!   do xp = 1,user%gridsize
!      do yp = 1,user%gridsize
!         call regolith_push(surf(xp,yp),cleanrego)
!      end do
!   end do
      
!   cleanrego%thickness = 10.0_DP !5000.0_DP
!   cleanrego%meltfrac = 0._DP
!   cleanrego%comp = 0._DP
!   do yp = 1,user%gridsize
!      do xp = 1,user%gridsize
!         call regolith_push(surf(xp,yp),cleanrego)
!      end do
!   end do

   !type(regolayertype) :: cleanrego
   !integer(I4B),parameter :: nlayers = 200
!   cleanrego%thickness = 20.0_DP !5000.0_DP
!   cleanrego%meltfrac = 0._DP
!   cleanrego%comp = 0._DP
!   do yp = 1,user%gridsize
!      do xp = 1,user%gridsize
!         do i=1,nlayers
!            call regolith_push(surf(xp,yp),cleanrego)
!         end do
!      end do
!   end do

! Mare/highlands contact 
   do yp = 1,user%gridsize
      do xp = 1,user%gridsize
         allocate(surf(xp,yp)%regolayer)
         nullify(surf(xp,yp)%regolayer%next)
      end do
   end do
 
   do yp = 1,user%gridsize
      do xp = 1,user%gridsize
!!         surf(xp,yp)%nmixi = 0
!!        surf(xp,yp)%nmixf = 0
!!         surf(xp,yp)%dexcav = 0._DP
!!         surf(xp,yp)%dmix = 0._DP
         if (xp<=user%gridsize/2) then 
            highland%thickness = 1000.0_DP
            highland%meltfrac  = 0._DP 
            highland%comp      = 0._DP
            call regolith_push(surf(xp,yp),highland)
            mare%thickness = 4000.0_DP
            mare%meltfrac  = 0._DP
            mare%comp      = 1.0_DP
            call regolith_push(surf(xp,yp),mare)
         else  
            highland%thickness = 5000.0_DP
            highland%meltfrac  = 0._DP
            highland%comp      = 0.0_DP
            call regolith_push(surf(xp,yp),highland)
         end if
      end do
   end do

! Square in the center
!  maresize = 300
!  do yp = 1,user%gridsize
!     do xp = 1,user%gridsize
!        if (xp>=user%gridsize/2-maresize .and. xp<=user%gridsize/2+maresize .and. &
!            yp>=user%gridsize/2-maresize .and. yp<=user%gridsize/2+maresize) then
!            highland%thickness = 4000.0_DP
!            highland%meltfrac  = 0._DP 
!            highland%comp      = 0._DP
!            call regolith_push(surf(xp,yp),highland)
!            mare%thickness = 1000.0_DP
!            mare%meltfrac  = 0._DP
!            mare%comp      = 1.0_DP
!            call regolith_push(surf(xp,yp),mare)
!        else
!            highland%thickness = 5000.0_DP
!            highland%meltfrac  = 0._DP
!            highland%comp      = 0._DP
!            call regolith_push(surf(xp,yp),highland)
!        end if
!     end do
!  end do

! Testing subcrater diffusion or subpixel regolith accumulation
!!   do yp=1,user%gridsize
!!      do xp=1,user%gridsize
!!            surf(xp,yp)%nmix = 0
            !surf(xp,yp)%dexcav = 0._DP
            !surf(xp,yp)%dmix = 0._DP
!!            highland%thickness = 200.0_DP
!!            highland%meltfrac  = 0.0_DP 
!!            highland%comp      = 0.0_DP
!!            call regolith_push(surf(xp,yp),highland)
            !mare%thickness = 40.0_DP
            !mare%meltfrac  = 0.0_DP
            !mare%comp      = 1.0_DP
            !call regolith_push(surf(xp,yp),mare)
!!      end do
!!   end do

! Testing subpixel mare  and highland contact
!   do yp=1,user%gridsize
!      do xp=1,user%gridsize
!         surf(xp,yp)%nmix = 0
!         surf(xp,yp)%dexcav = 0._DP
!         surf(xp,yp)%dmix = 0._DP
!         if (xp<=user%gridsize/2) then
!            highland%thickness = 2000.0_DP
!            highland%meltfrac  = 0.0_DP
!            highland%comp      = 0.0_DP
!            call regolith_push(surf(xp,yp),highland)
!            mare%thickness = 5.0_DP
!            mare%meltfrac  = 0.0_DP
!            mare%comp      = 1.0_DP
!            call regolith_push(surf(xp,yp),mare)
!         else
!            highland%thickness = 2005.0_DP
!            highland%meltfrac  = 0.0_DP
!            highland%comp      = 0.0_DP
!            call regolith_push(surf(xp,yp),highland)
!         end if 
!      end do
!   end do


   return
end subroutine init_regolith_stack

