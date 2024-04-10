!**********************************************************************************************************************************
!
!  Unit Name   : module_thermal
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Parameters and subroutine interface blocks for thermal routines 
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_thermal
use module_globals
implicit none
public
save

    interface
        subroutine thermal_initial_temperature(user,crater,thermi,distance)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(cratertype),intent(in) :: crater
        type(thermaltype),intent(inout) :: thermi
        real(DP),intent(in) :: distance
        end subroutine thermal_initial_temperature
    end interface

    interface
        subroutine thermal_dist(user,thermal,crater)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(cratertype),intent(in) :: crater
        end subroutine thermal_dist
    end interface

    interface
        subroutine thermal_diffusion(user,thermal)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        end subroutine thermal_diffusion
    end interface

end module