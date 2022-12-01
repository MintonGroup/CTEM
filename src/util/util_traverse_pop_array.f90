!**********************************************************************************************************************************
!
!  Unit Name   : util_traverse_pop_array
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
subroutine util_traverse_pop_array(regolayer,traverse_depth,poppedarray)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_traverse_pop_array
    implicit none
 
    ! Arguments
    !type(regolisttype),pointer   :: regolayer
    type(regodatatype),dimension(:),allocatable,intent(inout) :: regolayer
    real(DP),intent(in)          :: traverse_depth
    !type(regolisttype),pointer   :: poppedlist 
    type(regodatatype),dimension(:),allocatable,intent(out) :: poppedarray
 
    ! Internal variables
    real(DP)                    :: z,depth,dz
    type(regodatatype)          :: oldregodata
    logical                     :: initstat
    real(DP)                    :: recyratio
    integer(I4B)                :: i, N, maxi
 
    N = size(regolayer)
    depth = regolayer(N)%thickness
    dz = 0._DP
    z = traverse_depth
    
 
  
    ! Initialize popped array
    call util_init_array(poppedarray,initstat)
 
    
    !if (initstat) then
    i = N
    do
        depth = depth + regolayer(i)%depth
        if (depth > traverse_depth) then
            maxi = i
            exit
        else
            i = i - 1
        end if
    end do


    !allocate(poppedarray,source=regolayer(maxi:N))


    depth = regolayer(maxi)%depth

    ! if (z <= depth) then
    !     dz = depth - z
    oldregodata                  = regolayer(maxi)
    oldregodata%thickness        = z
    oldregodata%age(:)           = z / regolayer(maxi)%thickness * regolayer(maxi)%age(:)
    recyratio                    = dz / regolayer(maxi)%thickness
    regolayer(maxi)%age(:)    = recyratio * regolayer(maxi)%age(:)
    regolayer(maxi)%thickness = dz
    call util_push_array(poppedarray,oldregodata)
    ! else
    !     z = z - regolayer%thickness
    !     call util_pop_array(regolayer,oldregodata)
    !     call util_push_array(poppedarray,oldregodata)
    !     depth = regolayer%thickness
    ! end if
       
    return
 end subroutine util_traverse_pop_array
 