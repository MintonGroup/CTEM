!**********************************************************************************************************************************
!
!  Unit Name   : init_surf
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initializes the 3D thermal array, including background gradient
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine init_surf(user,thermal)
    use module_globals
    use module_init, EXCEPT_THIS_ONE => init_surf
    implicit none
 
    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:),intent(out) :: thermal

    ! Internal variables
    integer(I4B) :: i,j,k


    ! Executable code
    do i=1,user%gridsize
        do j=1,user%gridsize
            do k=1,user%zgridsize
                !calculate pixel depth
                thermal(i,j,k)%depth = (k-1) * user%zpix !this is depth at top of grid; can modify later
                !Add geothermal gradient
                thermal(i,j,k)%background = (13._DP/1000._DP) * thermal(i,j,k)%depth !13K/km for now; eventually could make it change over time?
                thermal(i,j,k)%temperature = thermal(i,j,k)%background
            end do
        end do
    end do
    return
end subroutine init_surf

