!****h* util/module_util
! Name
!   module_util -- Module for utility subroutines
! SYNOPSIS
!   This uses
!    NONE
!
!   use module_util
!
! DESCRIPTION
!
!   Miscellaneous utilities
!
! NOTES
!
!***

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

! interface
!    subroutine util_push(regolayer,newregodata)
!    use module_globals
!    implicit none
!    type(regolisttype),pointer :: regolayer
!    type(regodatatype),intent(in) :: newregodata
!    end subroutine util_push
! end interface

interface
   subroutine util_push_array(regolayer,newregodata)
   use module_globals
   implicit none
   type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
   type(regodatatype),intent(in) :: newregodata
   end subroutine util_push_array
end interface

interface
   subroutine util_pop_array(regolayer,oldregodata)
   use module_globals
   implicit none
   type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
   type(regodatatype),intent(out) :: oldregodata
   end subroutine util_pop_array
end interface


interface
   subroutine util_pop(regolayer,oldregodata)
   use module_globals
   implicit none
   type(regolisttype),pointer :: regolayer
   type(regodatatype),intent(out) :: oldregodata
   end subroutine util_pop
end interface

! interface
!    subroutine util_traverse_pop(regolayer,traverse_depth,poppedlist)
!    use module_globals
!    implicit none
!    !type(regolisttype),pointer :: regolayer
!    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
!    real(DP),intent(in)         :: traverse_depth
!    !type(regolisttype),pointer :: poppedlist
!    !type(regodatatype),dimension(:),allocatable,intent(out) :: poppedarray
!    end subroutine 
! end interface

interface
   subroutine util_traverse_pop_array(user,regolayer,traverse_depth,poppedarray)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
   real(DP),intent(in) :: traverse_depth
   type(regodatatype),dimension(:),allocatable,intent(out) :: poppedarray
   end subroutine
end interface

interface
   subroutine util_destroy_list(regolayer)
   use module_globals
   implicit none
   type(regolisttype),pointer :: regolayer
   end subroutine util_destroy_list
end interface

! interface
!    subroutine util_init_list(regolayer,initstat)
!    use module_globals
!    implicit none
!    type(regolisttype),pointer :: regolayer
!    logical, intent(out)     :: initstat
!    end subroutine util_init_list
! end interface

interface
   pure subroutine util_init_array(user,regolayer,domain,initstat)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
   type(domaintype),intent(in)    :: domain
   logical, intent(out)     :: initstat
   end subroutine util_init_array
end interface

interface
   subroutine util_add_to_layer(user,surfi,crater)
   use module_globals
   implicit none
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
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
   type(surftype),dimension(:,:),intent(inout) :: surf
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
   subroutine util_search_double(arr,ind,n,val,klo)
   use module_globals
   implicit none
   integer(I4B),intent(in) :: ind,n
   real(DP),dimension(:,:),intent(in) :: arr
   real(DP),intent(in) :: val
   integer(I4B),intent(inout) :: klo
   end subroutine util_search_double
   
   subroutine util_search_double_1(arr,ind,n,val,klo)
   use module_globals
   implicit none
   integer(I4B),intent(in) :: ind,n
   real(DP),dimension(:),intent(in) :: arr
   real(DP),intent(in) :: val
   integer(I4B),intent(inout) :: klo
   end subroutine util_search_double_1

   subroutine util_search_int(arr,ind,n,val,klo)
   use module_globals
   implicit none
   integer(I4B),intent(in) :: ind,n
   integer(I4B),dimension(:,:),intent(in) :: arr
   integer(I4B),intent(in) :: val
   integer(I4B),intent(inout) :: klo
   end subroutine util_search_int

end interface util_search

interface util_periodic
   pure subroutine util_periodic(x,y,side)
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
   subroutine util_push_regotemp(regotemp,regotime)
   use module_globals
   implicit none
   real(SP),dimension(:),allocatable,intent(inout) :: regotemp
   real(SP),dimension(:),allocatable,intent(inout) :: regotime
   end subroutine util_push_regotemp
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

interface
   function util_perlin_noise(x,y,z) result (noise)
   use module_globals
   implicit none
   real(DP),intent(in) :: x,y
   real(DP),intent(in),optional :: z
   real(DP) :: noise
   end function util_perlin_noise
end interface



! added by jundu on 10/25/2022
! generate random number with a normal distribution

interface
   subroutine util_random_number_uniform(u)
      use module_globals
      implicit none
      real(DP),intent(out) :: u
   end subroutine util_random_number_uniform
   subroutine util_random_number_normal(x)
      use module_globals
      implicit none
      real(DP),intent(out) :: x
   end subroutine util_random_number_normal
end interface











interface
   function util_npf_timefunc(T) result(N1)
   use module_globals
   real(DP), intent(in) :: T
   real(DP) :: N1
   end function util_npf_timefunc
end interface

interface
   function util_tscale(t) result(tscale)
   use module_globals
   real(DP), intent(in) :: t
   real(DP) :: tscale
   end function util_tscale
end interface

interface
   function util_t_from_scale(scale,start,finish) result(time)
   use module_globals
   real(DP), intent(in) :: scale, start, finish
   real(DP) :: time
   end function util_t_from_scale
end interface


end module

