!**********************************************************************************************************************************
!
!  Unit Name   : regolith_subcrater_mix
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
subroutine regolith_subcrater_mix(user,surf,domain,nflux,finterval,p)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_subcrater_mix
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(in) :: nflux ! impact rate (number of craters per m^2 per year)
   real(DP),intent(in) :: finterval  ! time elapsed ratio to the total time 
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

   ! Find the deepest depth for 100% true saturation
   if (p(2,1) < 1.0) then
      write(*,*) 'Error: The smallest sub-pixel craters has not covered one pixel-sized surface!'
   else
      klo = 1
      ds = p(1,1) 
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
                    dd = p(1,domain%smallest_impactor_index)
            end if

            if (surf(i,j)%regolayer%thickness < dd) then             
               call regolith_mix(surf(i,j),dd)
            end if
            
          end do
      end do 

   end if
   return
end subroutine regolith_subcrater_mix
