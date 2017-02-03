!**********************************************************************************************************************************
!
!  Unit Name   : util_search
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Implements hunt and bisection method search to find a position in an array for later interpolation
!  
!
!  Input
!    Arguments : arr  : 2D array ordered by dimension ind
!              : ind  : Index for ordered dimension
!              : n    : Size of array arr
!              : val  : Value to search such that val is between arr(k) and arr(k+1)
!              : k    : initial guess for the position
!
!  Output
!    Arguments : k    : one of three values: 1) position of val in arr. 
!                                            2) 0 = val is off the scale to the left
!                                            3) n+1 = val is off the scale to the right
!           
! 
!  Notes       :  Based on Numerical Recipes in Fortran 77 hunt function, with minor modifications
!
!**********************************************************************************************************************************
subroutine util_search_double(arr,ind,n,val,klo)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_search_double
implicit none

   ! Arguments
   integer(I4B),intent(in) :: ind,n
   real(DP),dimension(:,:),intent(in) :: arr
   real(DP),intent(in) :: val
   integer(I4B),intent(inout) :: klo

   ! Internals
   integer(I4B) :: k,khi,inc
   logical :: uporder ! true if arr is in ascending order

   ! Executable code
   uporder = arr(ind,1) <= arr(ind,n) 
   ! First hunt for the starting value
   if (.not.((klo <= 0).or.(klo > n))) then  ! Only proceed if guess is useful
      inc = 1
      if ((val >= arr(ind,klo)) .eqv. uporder) then
         do 
            khi = klo + inc 
            if (khi > n) then
               khi = n + 1
            else if ((val >= arr(ind,khi)) .eqv. uporder) then
               klo = khi
               inc = 2 * inc
               cycle
            end if
            exit
         end do
      else
         khi = klo
         do
            klo = khi - inc
            if (klo < 1 ) then
               klo = 0
            else if ((val < arr(ind,klo)) .eqv. uporder) then
               khi = klo
               inc = 2 * inc
               cycle
            end if
            exit
         end do
       end if
   else
      klo = 0 
      khi = n + 1
   end if

   ! Start the bisection method search
   do
      if (khi - klo == 1) then
         if (val == arr(ind,n)) klo = n-1
         if (val == arr(ind,1)) klo = 1
         return
      end if
      k = (khi + klo) / 2
      if ((val >= arr(ind,k)) .eqv. uporder) then
         klo = k
      else
         khi = k
      end if
   end do

end subroutine util_search_double

subroutine util_search_double_1(arr,ind,n,val,klo)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_search_double
implicit none

   ! Arguments
   integer(I4B),intent(in) :: ind,n
   real(DP),dimension(:),intent(in) :: arr
   real(DP),intent(in) :: val
   integer(I4B),intent(inout) :: klo

   ! Internals
   integer(I4B) :: k,khi,inc
   logical :: uporder ! true if arr is in ascending order

   ! Executable code
   uporder = arr(1) <= arr(n) 
   ! First hunt for the starting value
   if (.not.((klo <= 0).or.(klo > n))) then  ! Only proceed if guess is useful
      inc = 1
      if ((val >= arr(klo)) .eqv. uporder) then
         do 
            khi = klo + inc 
            if (khi > n) then
               khi = n + 1
            else if ((val >= arr(khi)) .eqv. uporder) then
               klo = khi
               inc = 2 * inc
               cycle
            end if
            exit
         end do
      else
         khi = klo
         do
            klo = khi - inc
            if (klo < 1 ) then
               klo = 0
            else if ((val < arr(klo)) .eqv. uporder) then
               khi = klo
               inc = 2 * inc
               cycle
            end if
            exit
         end do
       end if
   else
      klo = 0 
      khi = n + 1
   end if

   ! Start the bisection method search
   do
      if (khi - klo == 1) then
         if (val == arr(n)) klo = n-1
         if (val == arr(1)) klo = 1
         return
      end if
      k = (khi + klo) / 2
      if ((val >= arr(k)) .eqv. uporder) then
         klo = k
      else
         khi = k
      end if
   end do

end subroutine util_search_double_1


subroutine util_search_int(arr,ind,n,val,klo)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_search_int
implicit none

   ! Arguments
   integer(I4B),intent(in) :: ind,n
   integer(I4B),dimension(:,:),intent(in) :: arr
   integer(I4B),intent(in) :: val
   integer(I4B),intent(inout) :: klo

   ! Internals
   integer(I4B) :: k,khi,inc
   logical :: uporder ! true if arr is in ascending order

   ! Executable code
   uporder = arr(ind,1) <= arr(ind,n) 
   ! First hunt for the starting value
   if (.not.((klo <= 0).or.(klo > n))) then  ! Only proceed if guess is useful
      inc = 1
      if ((val >= arr(ind,klo)) .eqv. uporder) then
         do 
            khi = klo + inc 
            if (khi > n) then
               khi = n + 1
            else if ((val >= arr(ind,khi)) .eqv. uporder) then
               klo = khi
               inc = 2*inc
               cycle
            end if
            exit
         end do
      else
         khi = klo
         do
            klo = khi - inc
            if (klo < 1 ) then
               klo = 0
            else if ((val < arr(ind,klo)) .eqv. uporder) then
               khi = klo
               inc = 2*inc
               cycle
            end if
            exit
         end do
       end if
   else
      klo = 0 
      khi = n + 1
   end if
   ! Start the bisection method search

   do
      if (khi - klo == 1) then
         if (val == arr(ind,n)) klo = n-1
         if (val == arr(ind,1)) klo = 1
         return
      end if
      k = (khi+klo)/2
      if ((val >= arr(ind,k)) .eqv. uporder) then
         klo = k
      else
         khi = k
      end if
   end do

end subroutine util_search_int
