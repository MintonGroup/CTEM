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
subroutine thermal_dist(user,surf,thermal,crater)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_dist
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(surftype),dimension(:,:),intent(in) :: surf
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater

    ! Internal variables
    integer(I4B) :: inc,i,j,xpi,ypi,zinc,k,n
    real(DP) :: iradsq,lradsq,xp,yp,x_relative,y_relative,distance,lradcubed

    ! Executable code

    inc = max(min((13*crater%imprad/user%pix),real(user%gridsize-1,kind=DP)),1.0_DP) +1 ! setting this to 13*impactor radius for now; 
                                                                           ! calculations show the temperature increase is ~10 K.
    zinc = max(min((13*crater%imprad/user%zpix),real(user%zgridsize-1,kind=DP)),1.0_DP) +1

    do j=-inc,inc
        do i = -inc,inc
            iradsq = i**2 + j**2
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
            xp = xpi*user%pix
            yp = ypi*user%pix

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            x_relative = (crater%xl - xp)
            y_relative = (crater%yl - yp)
            lradsq = x_relative**2 + y_relative**2
            lradcubed = lradsq + surf(xpi,ypi)%dem

            n = 0
            
            do k=1,user%zgridsize
                if (n <= zinc) then
                    !calculate the 3-dimensional distance from layer depth
                    if (thermal(xpi,ypi,k)%depth > 0) then
                        distance = sqrt(lradsq+(thermal(xpi,ypi,k)%depth**2)) !currently uses "top left" instead of midpoint..
                        if (distance == 0.0_DP) then
                            distance = 750. !this could be avoided by using midpoint. 750 is just a test for now
                        end if
                        call thermal_initial_temperature(user,crater,thermal(xpi,ypi,k),distance)
                        n = n + 1
                    else
                        continue
                    end if
                else
                    exit
                end if
            end do
            
        end do
    end do

    return
end subroutine thermal_dist