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
        subroutine thermal_dist(user,surf,thermal,crater)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(surftype),dimension(:,:),intent(in) :: surf
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(cratertype),intent(in) :: crater
        end subroutine thermal_dist
    end interface

    interface
        subroutine thermal_diffusion(user,thermal,domain,difftime,icrater)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(domaintype),intent(in) :: domain
        real(DP), intent(in) :: difftime
        integer(I4B),intent(in) :: icrater
        end subroutine thermal_diffusion
    end interface

    interface
        subroutine thermal_depth_calculation(user,surf,crater,domain,thermal)
        use module_globals
        type(usertype),intent(in) :: user
        type(surftype),dimension(:,:),intent(in) :: surf
        type(cratertype),intent(in) :: crater
        type(domaintype),intent(inout) :: domain
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        end subroutine thermal_depth_calculation
    end interface

    interface
        subroutine thermal_ejecta_average(user,surf,crater,thermal,avgtemp)
        use module_globals
        type(usertype),intent(in) :: user
        type(surftype),dimension(:,:),intent(in) :: surf
        type(cratertype),intent(in) :: crater
        type(thermaltype),dimension(:,:,:),intent(in) :: thermal
        real(DP),intent(out) :: avgtemp
        end subroutine thermal_ejecta_average
    end interface

    interface
        subroutine thermal_add_ejecta(user,thermal,crater,avgtemp,thickness,tx,ty)
        use module_globals
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(cratertype),intent(in) :: crater
        real(DP),intent(in) :: avgtemp
        real(DP),intent(in) :: thickness
        integer(I4B) :: tx, ty
        end subroutine thermal_add_ejecta
    end interface

    interface
        subroutine thermal_remove(user,thermal,crater,prev)
        use module_globals
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(cratertype),intent(in) :: crater
        real(DP),dimension(:,:,:),allocatable,intent(in) :: prev
        end subroutine thermal_remove
    end interface

    interface
        subroutine thermal_interior(user,thermal,crater,incval,nmeltsheet,vmeltsheet)
        use module_globals 
        implicit none
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(cratertype),intent(in) :: crater
        integer(I4B),intent(in)     :: incval, nmeltsheet
        real(DP),intent(in)         :: vmeltsheet
        end subroutine thermal_interior
    end interface

    interface
        subroutine thermal_warp(user,thermal,crater)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
        type(cratertype),intent(in) :: crater
        end subroutine thermal_warp
    end interface

    interface
        subroutine thermal_link(user,crater,thermal,surf)
        use module_globals
        implicit none
        type(usertype),intent(in) :: user
        type(cratertype),intent(in) :: crater
        type(thermaltype),dimension(:,:,:),intent(in) :: thermal
        type(surftype),dimension(:,:),intent(inout) :: surf
        end subroutine thermal_link
    end interface
    
    ! interface
    !     subroutine thermal_uplift(user,thermal,crater,oldtemps)
    !     use module_globals
    !     type(usertype),intent(in) :: user
    !     type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    !     type(cratertype),intent(in) :: crater
    !     type(thermaltype),dimension(:,:,:), intent(in) :: oldtemps
    !     end subroutine thermal_uplift
    ! end interface
    
    
end module