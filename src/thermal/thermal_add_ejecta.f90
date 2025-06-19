!**********************************************************************************************************************************
!
!  Unit Name   : thermal_add_ejecta.f90
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Adds a layer of ejecta with a given average temperature to the thermal model.
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       : The uppermost elevation must be the baseline for the array 
!
!**********************************************************************************************************************************
subroutine thermal_add_ejecta(user,thermal,crater,avgtemp,thickness,tx,ty)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_add_ejecta
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(cratertype),intent(in) :: crater
    real(DP),intent(in) :: avgtemp
    real(DP),intent(in) :: thickness
    integer(I4B) :: tx, ty

    ! Internal variables
    real(DP) :: voxtemp, rem
    integer(I4B) :: npix, i, k

    ! Executable code

    if (thickness < user%zpix) then ! Ejecta layer less than z-pixel; will average
        voxtemp = ((avgtemp * (thickness/user%zpix)) + (thermal(tx,ty,1)%temperature*((user%zpix-thickness)/(user%zpix))) / 2.0_DP)
    end if
    npix = int(thickness / user%zpix)
    rem = mod(thickness,user%zpix)
    ! Shift thermal distribution down by the npix
    k = 0
    do i=1,user%zgridsize
        if (thermal(tx,ty,i)%depth > 0) then
            k = k + 1
            if (thickness < user%zpix) then
                thermal(tx,ty,i)%temperature = thermal(tx,ty,i)%temperature + voxtemp
                exit
            else
                if (k <= npix) then
                    if (k == npix) then
                        thermal(tx,ty,i)%temperature = thermal(tx,ty,i)%temperature + ((avgtemp * (rem/user%zpix)) + (thermal(tx,ty,i)%temperature*((user%zpix-thickness)/(user%zpix))) / 2.0_DP)
                    else
                        thermal(tx,ty,i)%temperature = thermal(tx,ty,i)%temperature + avgtemp
                    end if
                else
                    if (i+npix <= user%zgridsize) then
                        thermal(tx,ty,i)%temperature = thermal(tx,ty,i+npix)%temperature
                    else !temperature is equal to the background of the deepst voxel
                        thermal(tx,ty,i)%temperature = thermal(tx,ty,user%zgridsize)%background
                    end if
                end if
            end if
        end if
    end do
        

    return
end subroutine thermal_add_ejecta
        
