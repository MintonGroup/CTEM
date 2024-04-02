!**********************************************************************************************************************************
!
!  Unit Name   : regolith_depth_calculation.f90
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculate depth of each layer assuming the top of the layer is z=0
!  
!
!  Input
!    Arguments : regolayer :: array
!
!  Output
!    Arguments :
!           
! 
!  Notes       : 
!
!**********************************************************************************************************************************
subroutine regolith_depth_calculation(regolayer)
    use module_globals
    use module_regolith, EXCEPT_THIS_ONE => regolith_depth_calculation
    implicit none

    ! Arguments
    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer

    ! Internal variables
    integer(I4B) :: i
    real(DP) :: d

    ! Executable code
    d = 0.0_DP
    regolayer(size(regolayer))%depth = d
    do i = size(regolayer)-1,1,-1
        regolayer(i)%depth = d + regolayer(i+1)%thickness
        d = regolayer(i)%depth
    end do

    return
end subroutine regolith_depth_calculation


