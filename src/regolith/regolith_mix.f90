!**********************************************************************************************************************************
!
!  Unit Name   : regolith_mix
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Simulate a lunar regolith reworking zone by vertical mixing of several layers in our push-pop system        
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : surf : Surface expression matrix
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_mix(user,surf,domain,nflux,finterval)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(in) :: nflux ! impact rate (number of craters per m^2 per year)
   real(DP),intent(in) :: finterval  ! time elapsed ratio to the total time 

   ! Probability function builtup Internals
   integer(I4B) :: i, j
   real(DP),dimension(2,domain%smallest_impactor_index) :: p ! probability for different mixing depth at a given time
   real(DP) :: telapsed, r, a_crat, t ! calculating nflux 
   real(DP) :: rn ! random number 
   real(DP) :: ds, dd  ! mixing depth for shallow and deep
   integer(I4B) :: klo, khi

   ! Regolith layers mixing internals
   type(regolayertype),pointer :: current
   logical  :: MIX, DMIX
   real(DP) :: z, z0, zmare, ztot
   type(regolayertype) :: snewlayer, dnewlayer

   telapsed = user%interval * finterval 
   write(*,*) telapsed

   do i=1,domain%smallest_impactor_index
      r = nflux(1,i) / 2.0_DP
      a_crat = PI * r**2
      t = nflux(2,i) * a_crat * telapsed
      p(1,i) = r
      p(2,i) = 1.0_DP - exp(-1.0_DP * t)
   end do

   ! Find the deepest depth for 100% true saturation
   if (p(2,1) < 1.0) then
      write(*,*) 'Error: The smallest sub-pixel craters has not covered one pixel-sized surface!'
   else
      klo = 1
      ds = nflux(1,1)/2.0_DP
      do i=2,domain%smallest_impactor_index
         if (p(2,i)<1.0_DP) then
             exit
         else
            klo = i
            ds = p(1,i)
         end if
      end do

      do j=1,user%gridsize
         do i=1,user%gridsize

            ! Simulate a random mixing depth by a given time elapsed
            !call random_number(rn)
            if (rn>=p(2,1)) then  !random number is 100% true saturation!
               dd = ds 
            else if (rn<=p(2,domain%smallest_impactor_index)) then !random number reaches the deepest depth in the file
                    dd = nflux(1,domain%smallest_impactor_index) / 2.0_DP
            else 
               ! Find the depth for random probability in general
               call util_search_double(p,2,domain%smallest_impactor_index,rn,klo)
               dd = p(1,klo)
            end if

            ! Start to mixing stuff up: 
            ! Input variables: 1) ds (depth that reaches 100% true saturation) 
            !                  2) dd (some probability digging down deeper)
            ! There are several things that we should test before I mix layers.
            ! Becuase we may have a thick layer so that the mixing takes no effect.
            ! Then, if the first layer at a given pixel is thinner than the mixing depth, we will consider mixing.
            ! However, we have two mixing depth to do.
            ! Once we point to the first layer, we setup two flags to tell if we need to do either shallow or 
            ! deep or both mixing: MIX and DMIX. 
            ! From there, we also need to record the mixing ratio for both mixing layers, since the order of two layers
            ! is important. That is, the deeper layer should be top on the shallower layer so that the deeper layer
            ! should be pushed first, and the shallower layer comes second. For easy purpose to push, I will calculate 
            ! the mixing ratios of both layers first and do push-pop later. More importantly, before pushing, the unmixed
            ! layer should be popped first. As you can see, storing those layers will be easy to manupulate those layers 
            ! manually.  

            current => surf(i,j)%regolayer
            z = current%thickness
            MIX  = .false. 
            !DMIX = .false. 
 
            !if (dd > z .and. abs(dd-ds)>VSMALL) then
            ! We confirm that layers are thinner thah our shallow layer, ds. 
            if (z < ds .and. z > VSMALL) then
               zmare = 0._DP
               ztot  = 0._DP
               z0    = 0._DP
               do 
                if (.not. associated(current%next)) exit
                MIX = .true.
                if (z<=ds) then
                   ztot  = ztot  + current%thickness
                   zmare = zmare + current%thickness * current%comp
                   current => current%next
                   z0 = z
                   z  = z + current%thickness
                else
                   ztot  = ztot  + (ds - z0)
                   zmare = zmare + (ds - z0) * current%comp
                   z0    = ds
                   exit
                end if
               end do
               call regolith_traverse_pop(-1.0_DP * ztot, surf(i,j))
            else if (z >= ds) then
                    MIX   = .true.
                    ztot  = ds
                    zmare = ds * current%comp
                    z0    = ztot
                    call regolith_traverse_pop(-1.0_DP * ztot, surf(i,j))
            end if
 
            !if (MIX .eqv. .true.) then
            !   snewlayer%thickness = ztot
            !   snewlayer%comp      = zmare / ztot
            !   snewlayer%meltfrac  = 0._DP
            !   call regolith_traverse_pop(-1.0_DP*ztot, surf(i,j))
            !   call regolith_push(surf(i,j),snewlayer)
            !end if
 
               !if (MIX .eqv. .true. .and. snewlayer%thickness /= ds) then
               !write(*,*) i,j,snewlayer%thickness,ds
               !end if
                
               ! If the first layer is too think for shallow mixing,
               ! the first layer is unmixed, but the bottom part of the first 
               ! layer may undergo the mixing by chance, so here we use "z0"
               ! as marker that we track of layers prior to the current layer.
               ! If the previous layer is still thicker than the shallow mixing
               ! depth, that means the depth difference needs to add it into 
               ! deeper mixing process.  

               !if (z > dd) then
               !   DMIX = .true. 
               !   ztot  = (dd - ds)
               !   zmare = ztot * current%comp
               !else
               !   DMIX = .true.
               !   ztot  = (z-z0)
               !   zmare = ztot *current%comp
               !   z0    = z
               !   current => current%next
               !   z = z + current%thickness 
               !   do 
               !    if (.not. associated(current%next)) exit
               !    if (z<=dd) then
               !       ztot  = ztot  + current%thickness
               !       zmare = zmare + current%thickness * current%comp
               !       current => current%next
               !       z0 = z
               !       z = z + current%thickness
               !    else
               !       ztot  = ztot  + (dd - z0)
               !       zmare = zmare + (dd - z0) * current%comp
               !       exit
               !    end if
               !   end do
               !end if
                     
               !if (DMIX .eqv. .true. .and. (ztot + snewlayer%thickness) /= dd) then
               !write(*,*) 'x',ztot,snewlayer%thickness,ztot + snewlayer%thickness,dd
               !end if 
               ! Probably roundoff error, about the 15th digit. Presumably, it is good enough.
               !if (DMIX .eqv. .true.) then
               !   dnewlayer%thickness = ztot
               !   dnewlayer%comp      = zmare / ztot
               !   dnewlayer%meltfrac  = 0._DP
               !end if

             !end if
 
             !if (MIX .eqv. .true.) then
             !   if (DMIX .eqv. .true.) then
             !      call regolith_traverse_pop(-1.0_DP*snewlayer%thickness, surf(i,j))
             !      call regolith_traverse_pop(-1.0_DP*dnewlayer%thickness, surf(i,j))
             !      call regolith_push(surf(i,j),dnewlayer)
             !      call regolith_push(surf(i,j),snewlayer)
             !   else  
             !      call regolith_traverse_pop(-1.0_DP*snewlayer%thickness, surf(i,j))
             !      call regolith_push(surf(i,j),snewlayer)
             !   end if
             !end if
   
          end do
      end do 

   end if

   return
end subroutine regolith_mix
