!**********************************************************************************************************************************
!
!  Unit Name   : util_toupper
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 2003
!
!  Description : Convert string to uppercase
!
!  Input
!    Arguments : string : string to convert
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : string : converted string
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_toupper(string)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_toupper(string)

! Modules
     USE module_globals
     USE module_util, EXCEPT_THIS_ONE => util_toupper
     IMPLICIT NONE

! Arguments
     CHARACTER(*), INTENT(INOUT) :: string

! Internals
     INTEGER(I4B) :: i, length, index

! Executable code
     length = LEN(string)
     DO i = 1, length
          index = IACHAR(string(i:i))
          IF ((index >= LOWERCASE_BEGIN) .AND. (index <= LOWERCASE_END)) THEN
               index = index + UPPERCASE_OFFSET
               string(i:i) = ACHAR(index)
          END IF
     END DO

     RETURN

END SUBROUTINE util_toupper

