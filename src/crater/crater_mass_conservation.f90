!**********************************************************************************************************************************
!
!  Unit Name   : crater_mass_conservation
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Adjusts the elevation of the domain in order to conserve mass  between crater emplacements
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
subroutine crater_mass_conservation(user,surf,crater)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_mass_conservation
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(in)  :: crater

   ! Internal variables
   integer(I4B) :: i,j,inc,xpi,ypi,startinc,endinc,incsq,iradsq,ntot
   real(DP) :: tdem
   logical :: resetflag

   ! Executable code
   if (crater%maxinc >= user%gridsize / 2) then 
      startinc = -user%gridsize / 2 
      endinc = user%gridsize / 2 - 1
      resetflag = .true.
   else
      startinc = -crater%maxinc
      endinc = crater%maxinc
      resetflag = .false.
      incsq = endinc**2
   end if

   tdem = 0._DP
   ntot = 0
   do i = startinc, endinc
      do j = startinc, endinc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j
         call util_periodic(xpi, ypi, user%gridsize)
         if (resetflag) then
            tdem = tdem + surf(xpi, ypi)%dem 
            ntot = ntot + 1
         else
            iradsq = i**2 + j**2
            if (iradsq <= incsq) then 
               tdem = tdem + surf(xpi, ypi)%dem  - surf(xpi, ypi)%demOrig
               ntot = ntot + 1
            end if
         end if
      end do
   end do

   tdem = tdem / ntot

   ! Modify the surface elevation and ejcov         
   do i = startinc, endinc
      do j = startinc, endinc
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j
         call util_periodic(xpi, ypi, user%gridsize)
         iradsq = i**2 + j**2
         if ((.not.resetflag).and.(iradsq > incsq)) cycle
         surf(xpi, ypi)%dem = surf(xpi, ypi)%dem - tdem
         surf(xpi, ypi)%ejcov = surf(xpi, ypi)%ejcov - tdem
         surf(xpi, ypi)%demOrig = surf(xpi, ypi)%dem 
      end do
   end do
   
   return
end subroutine crater_mass_conservation

