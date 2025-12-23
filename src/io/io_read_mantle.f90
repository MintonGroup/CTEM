!**********************************************************************************************************************************
!
!  Unit Name   : io_read_mantle
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Reads in files for emplacing pre-existing mantle ejecta
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
subroutine io_read_mantle(user,surf,domain)
    use module_globals
    use module_io, EXCEPT_THIS_ONE => io_read_mantle
    implicit none

    ! Arguments

    type(usertype),intent(in)                :: user
    type(surftype),dimension(:,:),intent(in) :: surf
    type(domaintype),intent(in) :: domain

    ! Internals

    !real(DP) :: x,y,z

    ! Executable code

    !open(999, 'SPA_ejecta_thickness_angle30dimp260vimp13.txt')

    !read(999, *) x,y,z

    !print *, x,y,z



end subroutine io_read_mantle