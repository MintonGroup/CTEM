!**********************************************************************************************************************************
!
!  Unit Name   : thermal_loss_calc
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Solves the 3D heat equation using the Finite Difference Method for timescales relevant to fractional loss, then calculates fractional loss where relevant
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
!  Notes       : -Currently uses the whole grid
!
!**********************************************************************************************************************************
subroutine thermal_loss_calc(user,crater,thermal,surf)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_loss_calc
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(cratertype),intent(in) :: crater
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(surftype),dimension(:,:),intent(inout) :: surf

    ! Internal variables
    real(DP) :: kappa, gamma, delta_t, top, bottom, term1, term2, term3, ts, f, t, dr2
    integer(I4B) :: i,j,k,x,y,z,time,maxtime,nchanged,xplusone,xminusone,yplusone,yminusone
    type(thermaltype),dimension(:,:,:),allocatable :: initial, prev !Temperature at previous timestep (to prevent "new" temperature values from being used in diffusion)
    real(DP),dimension(:,:,:),allocatable :: times, losses

    allocate(prev,source=thermal)
    allocate(initial,source=thermal)
    allocate(losses(user%gridsize,user%gridsize,user%zgridsize))
    losses(:,:,:) = -1.0_DP
    allocate(times(user%gridsize,user%gridsize,user%zgridsize))
    times(:,:,:) = -1.0_DP

    maxtime = 10000
    delta_t = 1e3_DP
    kappa = 1e-6_DP

    top = 0.0_DP !Temperature at top of stack
    bottom = thermal(1,1,user%zgridsize)%background !For now make it equal to the geothermal gradient value at the bottom voxel
    !write(*,*) "Doing diffusion for", maxtime, "timesteps."

    do time = 1,maxtime
        nchanged = 0
        do k = 1,user%zgridsize
            do j = 1,user%gridsize
                do i = 1,user%gridsize
                    x = i
                    y = j
                    z = k
                    ! Factor in the repeating boundary conditions for the x and y dimensions
                    if (i == 1) then
                        xminusone = user%gridsize
                        xplusone = x + 1
                    else if (i == user%gridsize) then
                        xplusone = 1
                        xminusone = x - 1
                    else
                        xplusone = x + 1
                        xminusone = x - 1
                    end if

                    if (j == 1) then
                        yminusone = user%gridsize
                        yplusone = y + 1
                    else if (j == user%gridsize) then
                        yplusone = 1
                        yminusone = y - 1
                    else
                        yplusone = y + 1
                        yminusone = y - 1
                    end if

                    ! Actually do the diffusion, factoring in the boundary conditions in the z dimension
                    if (thermal(x,y,z)%temperature .eq. 0.0_DP .and. thermal(x,y,z)%background .eq. 0.0_DP) then !pixel is part of the boundary
                        continue
                    else
                        if (k == 1) then
                            term3 = ((prev(x,y,z+1)%temperature - (2*prev(x,y,z)%temperature) + top) / (user%zpix**2))
                        else if (k == user%zgridsize) then
                            term3 = ((bottom - (2*prev(x,y,z)%temperature) + prev(x,y,z-1)%temperature) / (user%zpix**2))
                        else
                            term3 = ((prev(x,y,z+1)%temperature - (2*prev(x,y,z)%temperature) + prev(x,y,z-1)%temperature) / (user%zpix**2))
                        end if
                        term1 = ((prev(xplusone,y,z)%temperature - (2*prev(x,y,z)%temperature) + prev(xminusone,y,z)%temperature) / (user%pix**2))
                        term2 = ((prev(x,yplusone,z)%temperature - (2*prev(x,y,z)%temperature) + prev(x,yminusone,z)%temperature) / (user%pix**2))
                        thermal(x,y,z)%temperature = prev(x,y,z)%temperature + delta_t * kappa * (term1 + term2 + term3)

                        if (thermal(x,y,z)%temperature .le. (thermal(x,y,z)%background+1.0_DP))  then
                            thermal(x,y,z)%temperature = thermal(x,y,z)%background
                        else
                            nchanged = nchanged+1
                        end if
                    end if
                    
                    if (initial(x,y,z)%temperature > 500.0_DP) then
                        if (losses(x,y,z) < 0.0_DP) then
                            times(x,y,z) = time * delta_t
                            if (thermal(x,y,z)%temperature < (10+thermal(x,y,z)%background) .or. thermal(x,y,z)%temperature < (2+prev(x,y,z)%temperature)) then
                                t = times(x,y,z)
                                dr2  = exp(-2.10_DP*(1e4_DP/thermal(x,y,z)%temperature)+8.05_DP)
                                f = ((6.0_DP/PI**(1.5_DP))*(((PI)**2.0_DP)*dr2*t)**(0.5_DP))-(((3/(PI**2.0_DP))*((PI)**2.0_DP)*dr2*t))
                                if (f < 0._DP .or. f > 0.85_DP) then
                                    f = 1.0_DP-(6.0_DP/PI**2.0_DP)*exp((-(PI)**2.0_DP)*dr2*t)
                                end if
                                if (f > 0._DP .and. f < 1._DP) then
                                    losses(x,y,z) = f
                                end if
                            end if
                        end if
                    end if
                end do
            end do
        end do

        prev(:,:,:)%temperature = thermal(:,:,:)%temperature

        if (nchanged == 0) exit !every voxel has cooled to the background temperature

    end do

    call thermal_loss_link(user,crater,thermal,surf,losses)
    deallocate(prev,initial,losses,times)

return
end subroutine thermal_loss_calc