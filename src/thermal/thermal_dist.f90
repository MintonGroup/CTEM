!**********************************************************************************************************************************
!
!  Unit Name   : thermal_dist
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Add the initial thermal distribution from an impact based on the scaling law found in Abramov et al. (2013)
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
subroutine thermal_dist(user,thermal,crater)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_dist
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater

    ! Internal variables
    integer(I4B) :: inc,i,j,xpi,ypi
    real(DP) :: iradsq,lradsq,xp,yp,x_relative,y_relative

    ! Executable code

    inc = max(min(13*crater%imprad,real(user%gridsize,kind=DP)),1.0_DP) +1 ! setting this to 13*impactor radius for now; 
                                                                           ! calculations show the temperature increase is ~10 K.
    do j=-inc,inc
        do i = -inc,inc
            iradsq = i**2 + j**2
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
            xp = xpi*user%pix
            yp = xpi*user%pix

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            x_relative = (crater%xl - xp)
            y_relative = (crater%yl - yp)
            lradsq = x_relative**2 + y_relative**2

            !call thermal_initial_temperature(user,crater,thermal,lradsq)
            !Decide how best to perform the Abramov calculations on a fixed grid: Its own subroutine or in this one? What about the data structure?
        end do
    end do

    return
end subroutine thermal_dist