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
subroutine util_init_array(user,regolayer,domain,initstat)
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
    ! if (.not. associated(regolayer)) then
    !    allocate(regolayer, STAT=allocstat)
    !    if (allocstat == 0) then
    !       initstat = .true.
    !       nullify(regolayer%next)
          ! regolayer%regodata%thickness = sqrt(VBIG) ! This generates a buffer layer that the model should never reach if the run is structured properly
          ! regolayer%regodata%comp = 0.0_DP
          ! regolayer%regodata%meltfrac = 0.0_DP
          ! regolayer%regodata%porosity = 0.0_DP
          ! regolayer%regodata%age(:)   = 0.0_SP
    if (allocated(regolayer)) deallocate(regolayer)
    allocate(regolayer(1),stat=allocstat)
    if (allocstat == 0) then
        initstat = .true.
        regolayer(1)%thickness = sqrt(VBIG) ! This generates a buffer layer that the model should never reach if the run is structured properly
        regolayer(1)%comp = 0.0_DP
        regolayer(1)%meltfrac = 0.0_DP
        regolayer(1)%porosity = 0.0_DP
        regolayer(1)%age(:)   = 0.0_SP
        regolayer(1)%ejm = 0.0_DP
        regolayer(1)%ejmf = 0.0_DP
        allocate(regolayer(1)%meltdist(domain%rcnum))
        regolayer(1)%meltdist(:) = 0.0_SP
        allocate(regolayer(1)%distvol(domain%rcnum))
        regolayer(1)%distvol(:) = 0.0_SP
        regolayer(1)%meltvolume = 0.0_DP
        regolayer(1)%totvolume = regolayer(1)%thickness * user%pix * user%pix
    end if
    !    else
    !       write(*,*) 'util_init_list: Initialization failed. Exhausted memory.'
    !    end if
    ! else
    !    write(*,*) 'util_init_list: Initialization failed. Regolayer already associated.'
    ! end if
 
    return
 end subroutine util_init_array