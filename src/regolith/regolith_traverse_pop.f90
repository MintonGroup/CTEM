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
subroutine regolith_traverse_pop(elchange,surfi,mixedregodata)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_traverse_pop
   implicit none

   ! Arguments
   real(DP),intent(in)            :: elchange
   type(surftype),intent(inout)   :: surfi
   type(regodatatype),intent(out) :: mixedregodata

   ! Internal variables
   real(DP)                    :: z,depth,dz
   type(regodatatype)          :: oldregodata

   !=======================================
   ! Get initial layer's info and the 
   ! desired info that we want to modify! 
   !=======================================
   depth = surfi%regolayer%regodata%thickness
   dz = 0._DP
   z = elchange
   mixedregodata%comp = 0.0_DP

   if (z < 0._DP) then

      do 
         if (.not. associated(surfi%regolayer)) then
            write(*,*) 'Major error in regolith_traverse_pop!'
            exit
         end if

         if (abs(z) <= depth) then
            dz = depth - abs(z)
            surfi%regolayer%regodata%thickness = dz
            mixedregodata%comp = mixedregodata%comp + dz * surfi%regolayer%regodata%comp
            mixedregodata%thickness = mixedregodata%thickness + dz
            exit
         else
            z = abs(z) - surfi%regolayer%regodata%thickness
            call util_pop(surfi,oldregodata)
            mixedregodata%comp = mixedregodata%comp + oldregodata%thickness * oldregodata%comp
            mixedregodata%thickness = mixedregodata%thickness + oldregodata%thickness
            depth = surfi%regolayer%regodata%thickness
         end if

      end do
      
   end if
   return
end subroutine regolith_traverse_pop
