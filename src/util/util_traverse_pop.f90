!****f* util/util_traverse_pop
! Name
!   util_traverse_pop -- Traversely pop off layers based on traverse depth (see DESCRIPTION). 
! SYNOPSIS
!   This uses 
!   * module_globals
!   * module_util
!   
!   call util_traverse_pop(regolayer,traverse_depth,poppedlist)
!
! DESCRIPTION
!    
!   In CTEM's linked list structure, a linked list has a physical depth. 
!   This subroutine takes a dpeth that will be popped off as long as the depth is above this depth.
!   Not only popping layers above a certain depth, but alaso reserving all popped layers. 
!   To do so, this subroutine will be:
!   * initializing a new linked list for popped layers,
!   * as long as the input layer's head is associated and above the traverse depth, popping a layer,
!   * and then pushing this popped layer to the popped list.
!
! ARGUMENTS
!   Input
!   * regolayer      -- pointer to the top of the regolith stack
!   * traverse_depth -- a depth relative the top of a layer
!   
!   Output
!   * regolayer      -- pointer to the top of the modified regolith stack 
!   * poppedlist     -- pointer to all popped layers.
! 
! NOTES
!   If the traverse depth is smaller than the current layer's thickness, the layer's thickness will 
!   be just modified. The other properties such as porosity or melt fraction should not be changed.
!   
!
!***

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
   logical                     :: initstat
   real(DP)                    :: recyratio

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
          oldregodata                  = regolayer%regodata
          oldregodata%thickness        = z
          oldregodata%age(:)           = z / regolayer%regodata%thickness * regolayer%regodata%age(:)
          recyratio                    = dz / regolayer%regodata%thickness
          regolayer%regodata%age(:)    = recyratio * regolayer%regodata%age(:)
          regolayer%regodata%thickness = dz
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
