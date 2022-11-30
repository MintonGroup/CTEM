!**********************************************************************************************************************************
!
!  Unit Name   : util_push_array
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pushes a new regolith block onto an old surface
!  
!
!  Input
!    Arguments : regolayer :: allocatable array
!                newregodata  :: new layer to push onto the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine util_push_array(regolayer,newregodata)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_push_array
    implicit none
 
    ! Arguments
    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
    type(regodatatype),intent(in) :: newregodata

    ! Internal variables
    type(regodatatype), dimension(:), allocatable :: newlayer
    integer(I4B) :: nold

    ! Executable code

    nold = size(regolayer)

    allocate(newlayer(nold+1))
    newlayer(1:nold) = regolayer(1:nold)
    newlayer(nold+1) = newregodata ! need to see if this adds to the top or bottom of the layer stack
    call move_alloc(newlayer, regolayer)
end subroutine util_push_array


