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
subroutine thermal_loss_calc(user,crater,thermal,surf,avgtemp,times,losses,timestamp)
    use module_globals
    use, intrinsic :: ieee_arithmetic
    use module_thermal, EXCEPT_THIS_ONE => thermal_loss_calc
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(cratertype),intent(in) :: crater
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(surftype),dimension(:,:),intent(inout) :: surf
    real(DP),intent(in) :: avgtemp
    real(DP),dimension(:,:,:),intent(inout) :: times,losses
    real(DP),intent(in) :: timestamp

    ! Internal variables
    real(DP) :: kappa, gamma, delta_t, top, bottom, term1, term2, term3, ts, f, t, dr2, diff_delta_t, mv, c, arg, e
    integer(I4B) :: i,j,k,x,y,z,time,maxtime,nchanged,xplusone,xminusone,yplusone,yminusone,  logdiff, diff
    type(thermaltype),dimension(:,:,:),allocatable :: initial, prev !Temperature at previous timestep (to prevent "new" temperature values from being used in diffusion)

    allocate(prev,source=thermal)
    allocate(initial,source=thermal)

    maxtime = 10000
    !maxtime = 50
    delta_t = 1e3_DP
    kappa = 1e-6_DP
    diff_delta_t = (1.0_DP/(2.0_DP * kappa)) * ((1.0_DP/(user%pix**2))+(1.0_DP/(user%pix**2))+(1.0_DP/(user%zpix**2)))**(-1.0_DP) 
    logdiff = int(log10(diff_delta_t))
    diff = 0


    !Find the closest time for the ejecta temperature to get fully reset

    if (avgtemp > 0) then
        do i = 1,logdiff+4
            t = 10**i
            dr2  = exp(-2.10_DP*(1e4_DP/(avgtemp+223))+8.05_DP) ! Assumes -50C for surface temperature; change to +273 for 0C
            f = ((6.0_DP/PI**(1.5_DP))*(((PI)**2.0_DP)*dr2*t)**(0.5_DP))-(((3/(PI**2.0_DP))*((PI)**2.0_DP)*dr2*t))
            if (f < 0._DP .or. f > 0.85_DP) then
                f = 1.0_DP-(6.0_DP/PI**2.0_DP)*exp((-(PI)**2.0_DP)*dr2*t)
            end if

            if (f == 1.0_DP) then
                delta_t = t / float(maxtime)
            else if (i == logdiff+4) then
                diff = 1
            else
                continue
            end if

        end do
    else
        diff = 0
    end if

    !TEST (remove when done)
    ! delta_t = 1e3_DP
    ! diff = 0
    !!!!!!!

    top = 0.0_DP !Temperature at top of stack
    bottom = thermal(1,1,user%zgridsize)%background !For now make it equal to the geothermal gradient value at the bottom voxel
    !write(*,*) "Doing diffusion for", maxtime, "timesteps."
    !write(*,*) "delta_t:", delta_t/(60*60*24*365), "yr. (", log10(delta_t), ")"
    mv = maxval(thermal(:,:,:)%temperature)

    if (diff == 0) then

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
                        
                        if (initial(x,y,z)%temperature > initial(x,y,z)%background) then
                            if (losses(x,y,z) < 0.1_DP) then
                                times(x,y,z) = time * delta_t
                                if (thermal(x,y,z)%temperature < (10+thermal(x,y,z)%background) .or. time == maxtime) then
                                    t = times(x,y,z)
                                    dr2  = exp(-2.10_DP*(1e4_DP/(initial(x,y,z)%temperature+223))+8.05_DP) ! Assumes -50C for surface temperature; change to +273 for 0C
                                    f = ((6.0_DP/PI**(1.5_DP))*(((PI)**2.0_DP)*dr2*t)**(0.5_DP))-(((3/(PI**2.0_DP))*((PI)**2.0_DP)*dr2*t))
                                    if (f < 0._DP .or. f > 0.85_DP) then
                                        if ( -((PI)**2.0_DP)*dr2*t < log(tiny(1.0_DP)) ) then
                                            f = 1.0_DP
                                        else
                                            f = 1.0_DP-(6.0_DP/PI**2.0_DP)*exp((-(PI)**2.0_DP)*dr2*t)
                                        end if
                                    end if
                                    if (f >= 0.1_DP .and. f <= 1._DP) then
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
    end if

    mv = maxval(losses(:,:,:))

    if (mv > 0.1) call thermal_loss_link(user,crater,thermal,surf,losses,timestamp)
    deallocate(prev,initial)!,losses,times)

return
end subroutine thermal_loss_calc