!**********************************************************************************************************************************
!
!  Unit Name   : util_pop
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Popped a old regolith block onto an old surface
!  
!
!  Input
!    Arguments : regolayer :: pointer to the top of the regolith stack
!                oldlayer  :: old layer to pop off of the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine util_pop_array(regolayer,oldregodata)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_pop_array
    implicit none

    ! Arguments
    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
    type(regodatatype),intent(out) :: oldregodata

    ! Internal variables
    type(regodatatype), dimension(:), allocatable :: newlayer
    integer(I4B) :: nold

    ! Executable code

    nold = size(regolayer)

    allocate(newlayer(nold-1),source=regolayer(1:nold-1))
    !newlayer(1:nold-1) = regolayer(1:nold-1) ! could also be 2:nold depending on if top or bottom is popped off
    oldregodata = regolayer(nold)
    call move_alloc(newlayer, regolayer)
end subroutine util_pop_array