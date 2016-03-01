!**********************************************************************************************************************************
!
!  Unit Name   : regolith_traverse_pop
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pop crater by traversing the linked list 
!  
!
!  Input
!    Arguments : cdepth :: a depth at a certain distance after emplacing a crater 
!
!  Output
!    Arguments : surf  
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_traverse_pop(elchange,surfi)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_traverse_pop
   implicit none

   ! Arguments
   real(DP),intent(in)         :: elchange
   type(surftype),intent(inout):: surfi

   ! Internal variables
   !type(regolayertype),pointer :: current
   real(DP)                    :: z,depth,dz

   ! Executable code

   depth = surfi%regolayer%thickness
   dz = 0._DP
   z = elchange

   if (z < 0._DP) then

      do 
       if (.not. associated(surfi%regolayer%next)) exit

       if (abs(z)<=depth) then
          dz = depth - abs(z)
          surfi%regolayer%thickness = dz
          exit
       else
          z = abs(z) - surfi%regolayer%thickness
          call regolith_pop(surfi)
          depth = surfi%regolayer%thickness
       end if

      end do
      
   end if

!   current => surfi%regolayer
!   depth = current%thickness
!   dz = 0._DP
!   z = elchange
!   if (z < 0._DP) then
!      do 
!       if (.not. associated(current%next)) exit
!
!       if (abs(z)<=depth) then
!          dz = depth - abs(z)
!          surfi%regolayer%thickness = dz
!          exit
!       else
!          z = abs(z) - surfi%regolayer%thickness
!          call regolith_pop(surfi)
!          depth = surfi%regolayer%thickness
!       end if
!
!      end do
!      
!   end if
      
   return
end subroutine regolith_traverse_pop
