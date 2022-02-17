!**********************************************************************************************************************************
!
!  Unit Name   : io_ejecta_table
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Outputs the ejecta table for a given crater
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine io_ejecta_table(crater,domain,ejb,ejtble,filename)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_ejecta_table
   implicit none

   ! Arguments
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   type(ejbtype),dimension(:),intent(in) :: ejb
   integer(I4B),intent(in) :: ejtble
   character(*),intent(in) :: filename

   ! Internal variables
   integer(I4B),parameter :: LUN=7
   integer(I4B) :: k

   ! Executable code
   open(LUN, FILE=filename, status='replace')
      write(LUN,'("# trad  = ",ES12.5, " frad = ",ES12.5)') crater%rad,crater%frad
      write(LUN,'("# ejrim = " ES12.5, " ejdis = ",ES12.5," imp = ",ES12.5)') crater%ejrim,crater%ejdis,crater%imp
      write(LUN,'(A63)') '# "r (m)"     "h (m)"      "v (m/s)"    "ang (deg)"  "erad (m)"'
      do k=1,ejtble 
         write(LUN,'(5(ES13.5E3,1X))') exp(ejb(k)%lrad),exp(ejb(k)%thick),sqrt(ejb(k)%vesq),ejb(k)%angle/DEG2RAD, &
                                     exp(ejb(k)%erad)
      end do
   close(LUN)


   return
end subroutine io_ejecta_table

