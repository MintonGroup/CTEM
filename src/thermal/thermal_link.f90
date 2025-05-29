!**********************************************************************************************************************************
!
!  Unit Name   : thermal_link
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Links temperature from thermal field to thermal history in surflayer
!  
!
!  Input
!    Arguments : 
!
!
!  Output
!    Arguments :
!           
! 
!  Notes       : 
!
!**********************************************************************************************************************************
subroutine thermal_link(user,thermal,surfi)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_link
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(in) :: thermal
    type(surftype),dimension(:,:),intent(inout) :: surfi

    ! Internal variables
    integer(I4B) :: i,j,k,regosize,test

    ! Executable code
    do i=1,user%gridsize
        do j=1,user%gridsize
            regosize = size(surfi(i,j)%regolayer)
            do k=1,regosize
                test = size(surfi(i,j)%regolayer(k)%thermalhist)
            end do
        end do
    end do
end subroutine thermal_link



