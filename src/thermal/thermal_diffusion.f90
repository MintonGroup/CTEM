!**********************************************************************************************************************************
!
!  Unit Name   : thermal_diffusion
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Solves the 3D heat equation using the Finite Difference Method
!  
!
!  Input
!    Arguments : thermal: thermal dataet
!
!
!  Output
!    Arguments : right now diffusion steps are outputted for movies; at release time they should not be written.
!           
! 
!  Notes       : -Currently uses the whole grid
!
!**********************************************************************************************************************************
subroutine thermal_diffusion(user,crater,thermal,surf,domain,difftime,icrater)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_diffusion
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(cratertype),intent(in) :: crater
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal
    type(domaintype),intent(in) :: domain
    type(surftype),dimension(:,:),intent(inout) :: surf
    real(DP),intent(in) :: difftime !in Ga
    integer(I4B),intent(in) :: icrater !for testing purposes only

    ! Internal variables
    real(DP) :: kappa, gamma, delta_t, top, bottom, term1, term2, term3, ts
    integer(I4B) :: i,j,k,x,y,z,time,maxtime,nchanged,xplusone,xminusone,yplusone,yminusone
    type(thermaltype),dimension(:,:,:),allocatable :: prev !Temperature at previous timestep (to prevent "new" temperature values from being used in diffusion)

    ! Test variables that will not be used in the actual release
    character(5) :: num
    character(19) :: filename

    ! Executable code

    kappa = 1e-6_DP !m/s^2; this is the value for "rock" (Jaeger et al., 1968; cited in Vaughn et al. 2013)
    delta_t = (1.0_DP/(2.0_DP * kappa)) * ((1.0_DP/(user%pix**2))+(1.0_DP/(user%pix**2))+(1.0_DP/(user%zpix**2)))**(-1.0_DP) !in s
    write(*,*) "delta_t:", delta_t/(60*60*24*365), "yr."

    if (user%testflag .eqv. .false. .or. domain%currentqmc .eqv. .true.) then
        ts = difftime * (60._DP * 60._DP * 24._DP * 365._DP * 1e9_DP)
        maxtime = ts / delta_t
    else
        maxtime = 10000 !diffusion test for testflag is an arbitrary number of timesteps
    end if

    !!!TEST DEBUG ONLY!!!!!!
    maxtime = 2
    !!!REMOVE THIS WHEN DONE!!!!!

    allocate(prev,source=thermal)

    !if(maxval(thermal(:,:,:)%temperature) > 500 .and. delta_t > 1e3) call thermal_loss_calc(user,crater,thermal,surf)

    top = 0.0_DP !Temperature at top of stack
    bottom = thermal(1,1,user%zgridsize)%background !For now make it equal to the geothermal gradient value at the bottom voxel
    write(*,*) "Doing diffusion for", maxtime, "timesteps."

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
                end do
            end do
        end do

        prev(:,:,:)%temperature = thermal(:,:,:)%temperature

        ! if (icrater == 1) then
        !     if (time == 1) then
        !         open(3,file='misc/therm00000.dat',status='replace',form='unformatted')
        !         write(3) prev(:,:,:)%temperature
        !         close(3)
        !     end if

        !     !prev(:,:,:)%temperature = thermal(:,:,:)%temperature

        !     ! Write out the timestep to the "misc" folder, which should be created already in the Python
        !     write(num,'(I0.5)') time
        !     filename = 'misc/therm'//trim(num)//'.dat'
        !     open(3,file=filename,status='replace',form='unformatted')
        !     write(3) thermal(:,:,:)%temperature
        !     close(3)
        ! else
        !     if (time == 1) then
        !         open(3,file='test/therm00000.dat',status='replace',form='unformatted')
        !         write(3) prev(:,:,:)%temperature
        !         close(3)
        !     end if

        !     !prev(:,:,:)%temperature = thermal(:,:,:)%temperature

        !     ! Write out the timestep to the "misc" folder, which should be created already in the Python
        !     write(num,'(I0.5)') time
        !     filename = 'test/therm'//trim(num)//'.dat'
        !     open(3,file=filename,status='replace',form='unformatted')
        !     write(3) thermal(:,:,:)%temperature
        !     close(3)
        ! end if

        if (nchanged == 0) exit !every voxel has cooled to the background temperature

    end do
    deallocate(prev)

return
end subroutine thermal_diffusion