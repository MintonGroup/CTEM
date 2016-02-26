!**********************************************************************************************************************************
!
!  Unit Name   : crater_make_list
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Makes a list of random indices from the production function to be used for generating the impactors 
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine crater_make_list(domain,prod,ntotcrat,production_list)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_make_list
   implicit none

   ! Arguments
   type(domaintype),intent(in) :: domain
   real(DP),dimension(:,:),intent(in) :: prod
   integer(I8B),intent(out) :: ntotcrat
   integer(I8B),dimension(:),intent(out) :: production_list 

   ! Internal variables
   integer(I4B) :: k
   real(DP) :: diffnum
 
   production_list = 0 
   do k = domain%smallest_impactor_index,domain%pnum
      if (k == domain%pnum) then
         diffnum = prod(2,k)
      else
         diffnum = prod(2,k) - prod(2,k+1)
      end if
      production_list(k) = util_poisson(diffnum)
   end do 
   ntotcrat = sum(production_list)

   return
end subroutine crater_make_list

