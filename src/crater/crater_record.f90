!**********************************************************************************************************************************
!
!  Unit Name   : crater_record
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : records the new crater in an available layer
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
subroutine crater_record(user,surf,crater)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_record
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater

   ! Internal variables
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   integer(I2B) :: isrim
   real(DP) :: rimdis
   

   ! Executable code

   ! determine area to effect
   rimdis = (crater%frad / user%pix)
   inc  = max(int(rimdis) + 1, 1)
   incsq = inc**2

   do j=-inc,inc
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! record crater in available layer
            call util_add_to_layer(user,surf(xpi,ypi),crater)
         end if

      end do
   end do 

   return
end subroutine crater_record

