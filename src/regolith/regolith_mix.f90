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
subroutine regolith_mix(user,surf,domain,nflux,p)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(in) :: nflux ! impact rate (number of craters per m^2 per year)
   real(DP),dimension(:,:),intent(in) :: p
   !real(DP),intent(out) :: ds

   ! Probability function builtup Internals
   integer(I4B) :: i, j
   real(DP) :: telapsed, r, a_crat, t ! calculating nflux 
   real(DP) :: rn ! random number 
   real(DP) :: dd, ds !mixing depth for shallow and deep
   integer(I4B) :: klo, khi

   ! Regolith layers mixing internals
   type(regolayertype),pointer :: current
   logical  :: MIX, DMIX
   real(DP) :: z, z0, zmare, ztot
   type(regolayertype) :: snewlayer, dnewlayer

   !do i=1,domain%smallest_impactor_index
   !   write(*,*) p(1,i), p(2,i)
   !end do

   ! Find the deepest depth for 100% true saturation
   if (p(2,1) < 1.0) then
      write(*,*) 'Error: The smallest sub-pixel craters has not covered one pixel-sized surface!'
   else
      klo = 1
      ds = nflux(1,1)/2.0_DP
      do i=2,domain%smallest_impactor_index
         if (p(2,i)<1.0_DP) exit
         klo = i
         ds = p(1,klo)
      end do

      do j=1,user%gridsize
         do i=1,user%gridsize

            call random_number(rn)
            if (rn < p(2,1) .and. rn > p(2,domain%smallest_impactor_index)) then
               call util_search_double(p,2,domain%smallest_impactor_index,rn,klo)
               dd = p(1,klo)
            else if (rn >= p(2,1)) then
                    dd = ds
            else
                    dd = nflux(1,domain%smallest_impactor_index) / 2.0_DP
            end if

            if (surf(i,j)%regolayer%thickness < dd) then 
            
               z = surf(i,j)%regolayer%thickness 
               zmare = 0._DP
               ztot  = 0._DP
               z0    = 0._DP
               do while ( (associated(surf(i,j)%regolayer%next)) .and. (z<dd) ) 
                     ztot  = ztot  + (z - z0)
                     zmare = zmare + (z - z0) * surf(i,j)%regolayer%comp
                     call regolith_pop(surf(i,j))
                     z0    = z
                     z     = z + surf(i,j)%regolayer%thickness
               end do

               ztot  = ztot  + (dd - z0)
               zmare = zmare + (dd - z0) * surf(i,j)%regolayer%comp
               call regolith_traverse_pop(-1.0_DP * (dd - z0), surf(i,j))
               dnewlayer%thickness = ztot 
               dnewlayer%comp      = zmare / ztot
               dnewlayer%meltfrac  = 0._DP
               call regolith_push(surf(i,j), dnewlayer)

            end if
            
          end do
      end do 

   end if
   return
end subroutine regolith_mix
