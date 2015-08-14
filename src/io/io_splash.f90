!**********************************************************************************************************************************
!
!  Unit Name   : io_splash
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes out splash text at the start of run
!
!  Input
!    Arguments : infile : input filename
!
!  Output
!    Arguments : 
! 
!  Notes       :  Subroutine sets several global variables based on user input file
!
!**********************************************************************************************************************************

subroutine io_splash()
   use module_globals
   !$ USE omp_lib
   use module_io, EXCEPT_THIS_ONE => io_splash
   implicit none

   ! Arguments

   ! Internals

   write(*,*) "*--------------------------------------------------------------------*"
   write(*,*) "*                       ____________________  ___                    *"
   write(*,*) "*                      / ____/_  __/ ____/  |/  /                    *"
   write(*,*) "*                     / /     / / / __/ / /|_/ /                     *"
   write(*,*) "*                    / /___  / / / /___/ /  / /                      *"
   write(*,*) "*                    \____/ /_/ /_____/_/  /_/                       *"
   write(*,*) "*                 Cratered Terrain Evolution Model                   *"
   write(*,*) "*--------------------------------------------------------------------*"
   write(*,*) " Version ",trim(adjustl(CTEMVER))
   write(*,*) " Authors: James E. Richardson <richardson@naic.edu>"
   write(*,*) "          David A. Minton <daminton@purdue.edu>" 
   write(*,*) "*--------------------------------------------------------------------*"

   !$ NTHREADS = OMP_get_max_threads() ! In the *parallel* case
   !$ write(*,'(a)')      ' OpenMP parameters:'
   !$ write(*,'(a)')      ' ------------------'
   !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads 
   !!$ write(*,*) 'Dynamic thread allocation: ',OMP_get_dynamic()
   return
end subroutine io_splash
