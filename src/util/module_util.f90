!**********************************************************************************************************************************
!
!  Unit Name   : util_toupper
!  Unit Type   : module
!  Project     : CTEM (adapted from Swifter by David E. Kaufmann)
!  Package     : util
!  Language    : Fortran 2003
!
!  Description : Miscellaneous utilities
!
!  Notes       : 
!
!**********************************************************************************************************************************
module module_util
implicit none
public
save

   interface
      subroutine util_push(surf,newlayer)
      use module_globals
      implicit none
      type(surftype),intent(inout) :: surf
      type(regodatatype),intent(in) :: newlayer
      end subroutine util_push
   end interface

   interface
      subroutine util_pop(surfi)
      use module_globals
      implicit none
      type(surftype),intent(inout):: surfi
      end subroutine util_pop
   end interface

interface
   subroutine util_add_to_layer(user,surfi,isrim,fcrat,xl,yl,depth,baseline)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   integer(I2B),intent(in) :: isrim
   real(DP),intent(in) :: fcrat
   real(SP),intent(in) :: xl,yl
   real(SP),intent(in) :: depth,baseline
   end subroutine util_add_to_layer
end interface

interface
   subroutine util_remove_from_layer(surfi,layer)
   use module_globals
   implicit none
   type(surftype),intent(inout) :: surfi
   integer(I4B),intent(in) :: layer
   end subroutine util_remove_from_layer
end interface


interface
   subroutine util_sort_layer(user,surf,crater)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   type(cratertype),intent(in) :: crater
   end subroutine util_sort_layer
end interface

interface
   SUBROUTINE util_toupper(string)
   USE module_globals
   IMPLICIT NONE
   CHARACTER(*), INTENT(INOUT) :: string
   END SUBROUTINE
end interface

interface
   Subroutine util_mrgrnk (XDONT, IRNGT)
   Use module_globals
   real(DP), dimension (:), intent (in) :: XDONT
   integer(I4B), dimension (:), intent (out) :: IRNGT
   end Subroutine util_mrgrnk
end interface

interface
   function util_rootbracketed(x1,y1) result(resultat)
   use module_globals
   implicit none
   real(DP),intent(in) :: x1,y1
   logical :: resultat
   end function util_rootbracketed
end interface


interface util_search
   subroutine util_search_double(arr,ind,n,val,k)
   use module_globals
   implicit none
   integer(I4B),intent(in) :: ind,n
   real(DP),dimension(:,:),intent(in) :: arr
   real(DP),intent(in) :: val
   integer(I4B),intent(out) :: k
   end subroutine util_search_double

   subroutine util_search_int(arr,ind,n,val,k)
   use module_globals
   implicit none
   integer(I4B),intent(in) :: ind,n
   integer(I4B),dimension(:,:),intent(in) :: arr
   integer(I4B),intent(in) :: val
   integer(I4B),intent(out) :: k
   end subroutine util_search_int

end interface util_search

interface util_periodic
   subroutine util_periodic(x,y,side)
   use module_globals
   implicit none
   integer(I4B),intent(inout) :: x,y
   integer(I4B),intent(in) :: side
   end subroutine util_periodic
end interface util_periodic

interface 
   function util_chi2test(data1,data2,nbins,p) result(ans)
   use module_globals
   implicit none
   integer(I4B),intent(in) :: nbins
   real(DP),dimension(nbins),intent(in) :: data1,data2
   real(DP),intent(out) :: p
   real(DP) :: ans
   end function util_chi2test
end interface

interface
   subroutine util_filecopy(file1,file2)
   use module_globals
   implicit none
   character(*),intent(in) :: file1,file2
   end subroutine util_filecopy
end interface

interface
   subroutine util_diffusion_solver(user,surf,N,indarray,kdiff,cumulative_elchange,maxhits)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   integer(I4B),intent(in) :: N
   real(DP),dimension(N,N),intent(in) :: kdiff
   integer(I4B),dimension(2,N,N),intent(in) :: indarray
   real(DP),dimension(N,N),intent(out) :: cumulative_elchange 
   integer(I4B),intent(in) :: maxhits
   end subroutine util_diffusion_solver
end interface

interface
   function util_area_intersection(R,xbar,ybar,P) result(area)
   use module_globals
   implicit none
   real(DP),intent(in) :: R,xbar,ybar,P
   real(DP) :: area
   end function util_area_intersection
end interface

interface
   function util_poisson(mu,poisson_first) result(ival)
   use module_globals
   implicit none
   real(DP), intent(in)    :: mu
   logical, intent(in),optional :: poisson_first
   integer(I8B)             :: ival
   end function util_poisson
end interface

end module

