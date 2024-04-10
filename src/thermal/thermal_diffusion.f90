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
subroutine thermal_diffusion(user,thermal)
    use module_globals
    use module_thermal, EXCEPT_THIS_ONE => thermal_diffusion
    implicit none

    ! Arguments
    type(usertype),intent(in) :: user
    type(thermaltype),dimension(:,:,:),intent(inout) :: thermal

    ! Internal variables
    real(DP) :: kappa, gamma, delta_t, top, bottom, term1, term2, term3
    integer(I4B) :: i,j,k,x,y,z,time,maxtime,nchanged,xplusone,xminusone,yplusone,yminusone
    type(thermaltype),dimension(:,:,:),allocatable :: prev !Temperature at previous timestep (to prevent "new" temperature values from being used in diffusion)

    ! Test variables that will not be used in the actual release
    character(5) :: num
    character(18) :: filename

    ! Executable code
    maxtime = 500 !number of timesteps; eventually may make this very high or perhaps devise a way to calculate. For now it's just a test

    kappa = 1e-6_DP !m/s^2; this is the value for "rock" (Jaeger et al., 1968; cited in Vaughn et al. 2013)
    delta_t = (1.0_DP/2.0_DP * kappa) * ((1.0_DP/(user%gridsize**2))+(1.0_DP/(user%gridsize**2))+(1.0_DP/(user%zgridsize**2)))**(-1.0_DP) !in s

    ! allocate(prev(user%gridsize,user%gridsize,user%zgridsize))
    ! prev(:,:,:)%temperature = thermal(:,:,:)%temperature
    allocate(prev,source=thermal)

    top = 0.0_DP !Temperature at top of stack
    bottom = thermal(1,1,1)%background !For now make it equal to the geothermal gradient value at the bottom voxel

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

                    ! Skip pixels where the temperature is already equal to the background
                    if (thermal(x,y,z)%temperature == thermal(x,y,z)%background) cycle

                    ! Actually do the diffusion, factoring in the boundary conditions in the z dimension
                    if (k == 1) then
                        term3 = ((prev(x,y,z+1)%temperature - (2*prev(x,y,z)%temperature) + bottom) / (user%zpix**2))
                    else if (k == user%zgridsize) then
                        term3 = ((top - (2*prev(x,y,z)%temperature) + prev(x,y,z-1)%temperature) / (user%zpix**2))
                    else
                        term3 = ((prev(x,y,z+1)%temperature - (2*prev(x,y,z)%temperature) + prev(x,y,z-1)%temperature) / (user%zpix**2))
                    end if
                    term1 = ((prev(xplusone,y,z)%temperature - (2*prev(x,y,z)%temperature) + prev(xminusone,y,z)%temperature) / (user%pix**2))
                    term2 = ((prev(x,yplusone,z)%temperature - (2*prev(x,y,z)%temperature) + prev(x,yminusone,z)%temperature) / (user%pix**2))
                    thermal(x,y,z)%temperature = prev(x,y,z)%temperature + delta_t * kappa * (term1 + term2 + term3)

                    if (abs(thermal(x,y,z)%temperature - thermal(x,y,z)%background) .lt. 1e-1) then !high tolerance for now; if needed make closer to 1e-16
                        thermal(x,y,z)%temperature = thermal(x,y,z)%background
                    else
                        nchanged = nchanged+1
                    end if
                    if (nchanged == 0) exit !every voxel has cooled to the background temperature
                end do
            end do
        end do
        prev(:,:,:)%temperature = thermal(:,:,:)%temperature

        ! Write out the timestep to the "misc" folder, which should be created already in the Python <--this is just for the test example
        write(num,'(I0.5)') time
        filename = 'misc/therm'//trim(num)//'.dat'
        open(3,file=filename,status='replace',form='unformatted')
        write(3) thermal(:,:,:)%temperature
        close(3)

    end do
    deallocate(prev)

return
end subroutine thermal_diffusion