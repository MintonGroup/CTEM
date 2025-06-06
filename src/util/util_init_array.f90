!**********************************************************************************************************************************
!
!  Unit Name   : util_init_array
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Initialize an new allocatable array
!  
!
!  Input
!    Arguments : regolayer :: array
!                oldlayer  :: old layer to pop off of the top of the stack
!
!  Output
!    Arguments :
!           
! 
!  Notes       : 
!
!**********************************************************************************************************************************
pure subroutine util_init_array(user,regolayer,domain,initstat)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_init_array
    implicit none
 
    ! Arguments
    type(usertype),intent(in) :: user
    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
    type(domaintype),intent(in)    :: domain
    logical, intent(out)       :: initstat
 
    ! Internal variables
    integer(I4B) :: allocstat 
 
    ! Executable code
    initstat = .false.
    if (allocated(regolayer)) deallocate(regolayer)
    allocate(regolayer(1),stat=allocstat)
    if (allocstat == 0) then
        initstat = .true.
        regolayer(1)%thickness = user%trad ! This generates a buffer layer that the model should never reach if the run is structured properly
        regolayer(1)%comp = 0.0_DP
        regolayer(1)%age(:)   = 0.0_SP
        regolayer(1)%ejm = 0.0_DP
        allocate(regolayer(1)%distvol(1+domain%rcnum))
        regolayer(1)%distvol(:) = 0.0_SP
        regolayer(1)%meltvolume = 0.0_DP
        regolayer(1)%totvolume = regolayer(1)%thickness * user%pix * user%pix
        allocate(regolayer(1)%regotemp(1,1))
        allocate(regolayer(1)%regotime(1,1))
        regolayer(1)%regotemp(1,1) = 0.0_SP
        regolayer(1)%regotime(1,1) = 0.0_SP
    end if
 
    return
 end subroutine util_init_array