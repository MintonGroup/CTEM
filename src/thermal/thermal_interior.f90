!**********************************************************************************************************************************
!
!  Unit Name   : thermal_interior.f90
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Adds an interior melt sheet to the threrml grid
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       : Assumes the melt sheet temperature is the coolest possible temperature of melt: 1087 C according to Abramov et al. (2013) 
!
!**********************************************************************************************************************************
subroutine thermal_interior(user,thermal,crater,incval,nmeltsheet,vmeltsheet)
    use module_globals 
    use module_util
    use module_thermal, EXCEPT_THIS_ONE => thermal_interior
    implicit none

    !Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater
    integer(I4B),intent(in)     :: incval, nmeltsheet
    real(DP),intent(in)         :: vmeltsheet

    !internal variables
    integer(I4B) xpi,ypi,i,j,k,inc,incsq,iradsq,npix
    real(DP) :: lradsq, x_relative, y_relative, xp, yp, hmeltsheet
    type(thermaltype),dimension(:,:,:),allocatable :: oldtemps

    !Executable code

    allocate(oldtemps,source=thermal(:,:,:))

    hmeltsheet = vmeltsheet / (nmeltsheet*user%pix*user%pix)
    npix = hmeltsheet / user%zpix

    inc = incval

    do j=-inc,inc
        do i=-inc,inc
            iradsq = i**2 + j**2

            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            xp = xpi * user%pix
            yp = ypi * user%pix

            call util_periodic(xpi,ypi,user%gridsize)
            x_relative = (crater%xl - xp)
            y_relative = (crater%yl - yp)
            lradsq = x_relative**2 + y_relative**2

            if (lradsq > crater%frad**2) cycle

            do k=1,user%zgridsize
                if (thermal(xpi,ypi,k)%depth > 0) then
                    if (hmeltsheet < user%zpix) then !Average melt sheet with the current temperature of the pixel
                        thermal(xpi,ypi,k)%temperature = max((1087._DP * (hmeltsheet / user%zpix))&
                         + (thermal(xpi,ypi,k)%temperature * ((user%zpix-hmeltsheet)/(user%zpix))),1087._DP)
                        exit
                    else !Emplace melt sheet of a given thickness, then push the rest of the pixels down
                        if (thermal(xpi,ypi,k)%depth < hmeltsheet) then
                            thermal(xpi,ypi,k)%temperature = 1087.
                        else
                            if (thermal(xpi,ypi,k)%depth - hmeltsheet < user%zpix) then !final pixel of melt sheet
                                thermal(xpi,ypi,k)%temperature = max((1087._DP * (thermal(xpi,ypi,k)%depth - hmeltsheet / user%zpix))&
                                    + (thermal(xpi,ypi,k)%temperature * ((user%zpix-(thermal(xpi,ypi,k)%depth - hmeltsheet))/(user%zpix))),1087._DP)
                            else !no melt sheet, but shift the distribution down
                                thermal(xpi,ypi,k)%temperature = oldtemps(xpi,ypi,k-npix)%temperature
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end do

    deallocate(oldtemps)

    return
    
end subroutine thermal_interior

