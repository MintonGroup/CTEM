!**********************************************************************************************************************************
!
!  Unit Name   : util_traverse_pop
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
!  Notes       :  Popped list will be in reversed order from the original list
!
!**********************************************************************************************************************************
subroutine util_traverse_pop(regolayer,traverse_depth,poppedlist)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_traverse_pop
   implicit none

   ! Arguments
   type(regolisttype),pointer   :: regolayer
   real(DP),intent(in)          :: traverse_depth
   type(regolisttype),pointer   :: poppedlist 

   ! Internal variables
   real(DP)                    :: z,depth,dz
   type(regodatatype)          :: oldregodata
   logical :: initstat

   depth = regolayer%regodata%thickness
   dz = 0._DP
   z = traverse_depth
   poppedlist => null()
 
   ! Initialize popped list
   call util_init_list(poppedlist,initstat)
   
   if (initstat) then
      do 
       if (.not. associated(regolayer)) then
          write(*,*) 'Major error in util_traverse_pop!'
          exit
       end if

       if (z <= depth) then
          dz = depth - z
          regolayer%regodata%thickness = dz
          oldregodata = regolayer%regodata
          oldregodata%thickness = z
          call util_push(poppedlist,oldregodata)
          exit
       else
          z = z - regolayer%regodata%thickness
          call util_pop(regolayer,oldregodata)
          call util_push(poppedlist,oldregodata)
          depth = regolayer%regodata%thickness
       end if
      
      end do
   else
      write(*,*) 'util_traverse_pop: Initialization of poppedlist failed.'
   end if

   return
end subroutine util_traverse_pop
