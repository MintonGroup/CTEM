!**********************************************************************************************************************************
!
!  Unit Name   : io_write_mantle
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in files for emplacing pre-existing mantle ejecta data
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine io_write_mantle(user,surf,domain)
    use module_globals
    use module_io, EXCEPT_THIS_ONE => io_write_mantle
    implicit none

    ! Arguments

    type(usertype),intent(in)                :: user
    type(surftype),dimension(:,:),intent(in) :: surf
    type(domaintype),intent(in) :: domain

    ! Internals

    character(len=*), parameter :: infile  = "SPA_ejecta_thickness_angle30dimp260vimp13.txt"
    character(len=*), parameter :: outdat = "mantle_ejecta.dat"

    real(8), parameter :: blank = 0.0d0
    !real(8), parameter :: xmin = -2000.0d0, xmax = 2000.0d0
    !real(8), parameter :: ymin = -2000.0d0, ymax = 2000.0d0

    integer :: iu_in, iu_out, ios
    character(len=1024) :: line
    real(8) :: x, y, c
    logical :: has_c

    ! Executable code

    open(newunit=iu_in, file=infile, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) "Error: cannot open ", infile
      stop 1
    end if

    open(newunit=iu_out, file=outdat, status="replace", action="write", iostat=ios)
    if (ios /= 0) then
      write(*,*) "Error: cannot open ", outdat
      stop 1
    end if

    do
      read(iu_in,'(A)', iostat=ios) line
      if (ios /= 0) exit

      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      ! Try 3 columns first
      read(line, *, iostat=ios) x, y, c
      if (ios == 0) then
        has_c = .true.
      else
        ! Try 2 columns
        read(line, *, iostat=ios) x, y
        if (ios /= 0) cycle
        has_c = .false.
        c = 1.0d0
      end if

      ! If there is a 3rd column and it is 0.0, treat as blank (skip)
      if (has_c) then
        if (c == blank) cycle
      end if

      ! Skip points outside the requested plotting window

      !if (x < xmin .or. x > xmax) cycle
      !if (y < ymin .or. y > ymax) cycle

      ! Write normal, mirrored, and thickness to one file

      write(iu_out,'(F16.6,1X,F16.6,1X,F16.6)') (x*1000/(user%pix))+((user%gridsize)/2), (y*1000/(user%pix))+((user%gridsize)/2), c*1000 ! converting km to pixels and shifting from (0,0) 
      write(iu_out,'(F16.6,1X,F16.6,1X,F16.6)') (x*1000/(user%pix))+((user%gridsize)/2), (-y*1000/(user%pix))+((user%gridsize)/2), c*1000 ! center to (gridsize/2, gridsize/2) center
      !write(iu_norm,'(F16.6,1X,F16.6,1X,F16.6)') y, x, c
      !write(iu_mirr,'(F16.6,1X,F16.6,1X,F16.6)') -y, x, c
    end do

    close(iu_in)
    close(iu_out)

    return
end subroutine io_write_mantle