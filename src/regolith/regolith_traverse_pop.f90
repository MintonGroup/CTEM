!**********************************************************************************************************************************
!
!  Unit Name   : regolith_traverse_pop
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Removes all layers down to a given depth. Cuts a layer if the depth ends in the middle of an old layer.
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : surf  
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine util_traverse_pop(elchange,surfi,poppedlist)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_traverse_pop
   implicit none

   ! Arguments
   real(DP),intent(in)            :: elchange
   type(surftype),intent(inout)   :: surfi
   type(regolisttype),pointer :: poppedlist => null()

   ! Internal variables
   real(DP)                    :: z,depth,dz
   type(regodatatype)          :: oldregodata

   !=======================================
   ! Get initial layer's info and the 
   ! desired info that we want to modify! 
   !=======================================
   depth = surfi%regolayer%thickness !we don't need '%regodata' anymore
   dz = 0._DP
   z = elchange
   mixedregodata%comp = 0.0_DP

   if (z < 0._DP) then

      do 
         ! if (.not. associated(surfi%regolayer)) then
         !    write(*,*) 'Major error in regolith_traverse_pop!'
         !    exit
         ! end if

         if (abs(z) <= depth) then
            dz = depth - abs(z)
            surfi%regolayer%thickness = dz
            mixedregodata%comp = mixedregodata%comp + dz * surfi%regolayer%comp
            mixedregodata%meltfrac = mixedregodata%meltfrac + dz * surfi%regolayer%%meltfrac
            mixedregodata%thickness = mixedregodata%thickness + dz
            exit
         else
            z = abs(z) - surfi%regolayer%regodata%thickness
            call util_pop_array(surfi,oldregodata)
            mixedregodata%comp = mixedregodata%comp + oldregodata%thickness * oldregodata%comp
            mixedregodata%meltfrac = mixedregodata%meltfrac + oldregodata%thickness * oldregodata%meltfrac
            mixedregodata%thickness = mixedregodata%thickness + oldregodata%thickness
            depth = surfi%regolayer%thickness
         end if

      end do
      
   end if
   return
end subroutine regolith_traverse_pop
