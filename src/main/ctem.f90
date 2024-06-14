!**********************************************************************************************************************************
!
!  Unit Name   : CTEM
!  Unit Type   : program
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : The Cratered Terrain Evolution Model
!
! 
!  Notes       :  This is the main function, which serves mainly to set up the
!  runs and handle input and output. Actually cratering is performed in the 
!  crater_populate subroutine
!
!**********************************************************************************************************************************
program CTEM
   use driver
   implicit none
   call ctem_driver()
end program
