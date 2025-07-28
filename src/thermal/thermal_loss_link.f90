!**********************************************************************************************************************************
!
!  Unit Name   : thermal_loss_link
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : links surf to fractional loss
!  
!
!  Input
!    Arguments : thermal: thermal dataet
!
!
!  Output
!    Arguments : 
!           
! 
!  Notes       : THE ARRAY "REGOTEMP" IS ACTUALLY FRACTIONAL LOSS!!!!   
!
!**********************************************************************************************************************************
subroutine thermal_loss_link(user,crater,thermal,surf,losses)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_loss_link
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(cratertype),intent(in) :: crater
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(surftype),dimension(:,:),intent(inout) :: surf
    real(DP),dimension(:,:,:),intent(in) :: losses

    integer(I4B) :: i,j,k,m,regosize,histsize

    real(DP) :: cum_thickness, current_depth, average_depth, current_therm_depth, tplusone, tminusone, temp
    real(DP) :: dplusone, dminusone, dpercent

    ! Executable code
    do j=1,user%gridsize
        do i=1,user%gridsize
            regosize = size(surf(i,j)%regolayer)
            cum_thickness = 0
            do k=1,regosize
                current_depth = cum_thickness + surf(i,j)%regolayer(k)%thickness 
                !interpolate between minimum (cum_thickness) and maximum (current_depth) depth for this layer
                average_depth = (cum_thickness + current_depth) / 2.0_DP
                do m=1,user%zgridsize !Find temperature at the thermal location corresponding to this depth
                    current_therm_depth = thermal(i,j,m)%depth
                    if (losses(i,j,m) > 0.0_DP) then
                        call util_push_regotemp(surf(i,j)%regolayer(k))
                        histsize = size(surf(i,j)%regolayer(k)%regotemp, 2)
                        surf(i,j)%regolayer(k)%regotemp(:,histsize) = losses(i,j,m)
                        surf(i,j)%regolayer(k)%regotime(:,histsize) = crater%timestamp
                    end if
                end do
            end do
        end do
    end do

end subroutine thermal_loss_link

