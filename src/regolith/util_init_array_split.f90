!**********************************************************************************************************************************
!
!  Unit Name   : util_init_array_split
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
pure subroutine util_init_array_split(user,regolayer,domain,initstat)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_init_array_split
implicit none

! Arguments
type(usertype),intent(in) :: user
type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
type(domaintype),intent(in)    :: domain
logical, intent(out)       :: initstat

! Internal variables
integer(I4B) :: allocstat, i
real(DP) :: therm_thickness

! Executable code
initstat = .false.
if (allocated(regolayer)) deallocate(regolayer)
allocate(regolayer(1+user%zgridsize),stat=allocstat)
therm_thickness = user%zgridsize * user%zpix
if (allocstat == 0) then
    initstat = .true.
    regolayer(1)%thickness = user%trad - therm_thickness ! This generates a buffer layer that the model should never reach if the run is structured properly
    regolayer(1)%comp = 0.0_DP
    regolayer(1)%age(:)   = 0.0_SP
    regolayer(1)%ejm = 0.0_DP
    allocate(regolayer(1)%distvol(1+domain%rcnum))
    regolayer(1)%distvol(:) = 0.0_SP
    regolayer(1)%meltvolume = 0.0_DP
    regolayer(1)%totvolume = regolayer(1)%thickness * user%pix * user%pix
    allocate(regolayer(1)%regotemp(1))
    allocate(regolayer(1)%regotime(1))
    regolayer(1)%regotemp(1) = 0.0_SP
    regolayer(1)%regotime(1) = 0.0_SP
    do i = 2,user%zgridsize+1
        regolayer(i)%thickness = user%zpix
        regolayer(i)%comp = 0.0_DP
        regolayer(i)%age(:)   = 0.0_SP
        regolayer(i)%ejm = 0.0_DP
        allocate(regolayer(i)%distvol(1+domain%rcnum))
        regolayer(i)%distvol(:) = 0.0_SP
        regolayer(i)%meltvolume = 0.0_DP
        regolayer(i)%totvolume = regolayer(1)%thickness * user%pix * user%pix
        allocate(regolayer(i)%regotemp(1))
        allocate(regolayer(i)%regotime(1))
        regolayer(i)%regotemp(1) = 0.0_SP
        regolayer(i)%regotime(1) = 0.0_SP
    end do
end if

return
end subroutine util_init_array_split