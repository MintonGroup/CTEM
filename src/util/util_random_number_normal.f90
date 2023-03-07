!****f* util/util_random_number_normal
! Name
!   util_random_number_normal -- Push a new layer onto the top of an old layer. 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   
!   call util_random_number_normal(x)
!
! DESCRIPTION
!    
!   Push subroutine is to create new head and push it with new data onto the old layer.
!   This subroutine will be: 
!   * checking if the head of an input old layer is associated.
!   * if we have enough space, allocating a new head and space for new data.  
!   * at the end, linking new head with new data to the old head. 
!
! ARGUMENTS
!   Input
!   * regolayer   -- pointer to the top of the regolith stack
!   * newregodata -- new regodata that is about to be pushed. 
!   
!   Output
!   * regolayer   -- pointer to the top of the regolith stack with newlayer. 
! 
! NOTES
! 
!***

!**********************************************************************************************************************************
!
!  Unit Name   : util_random_number_normal
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pushes a new regolith block onto an old surface
!  
!
!  Input
!    Arguments : regolayer :: pointer to the top of the regolith stack
!                newlayer  :: new layer to push onto the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************

subroutine util_random_number_uniform(u)
   use module_globals
   use module_util 
   implicit none
   real(DP),intent(out) :: u
   real(DP) :: r
   call random_number(r)
   u = 1 - r
end subroutine util_random_number_uniform

subroutine util_random_number_normal(x)
   use module_globals
   use module_util 
   implicit none
   real(DP),intent(out) :: x
   real(DP) :: u1,u2
   call util_random_number_uniform(u1)
   call util_random_number_uniform(u2)
   x = sqrt(-2*log(u1))*cos(2*PI*u2)
end subroutine util_random_number_normal