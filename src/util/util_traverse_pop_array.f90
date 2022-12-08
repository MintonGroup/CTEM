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
    type(regodatatype),dimension(:),allocatable         :: oldregodata
    !logical                     :: initstat
    real(DP)                    :: recyratio, depth_diff
    integer(I4B)                :: i, N, maxi
 
    N = size(regolayer)
    depth = 0._DP
    dz = 0._DP
    z = traverse_depth
    
    i = N
    ! do !this is where i=N,1,-1 could be used
    !     depth = depth + regolayer(i)%thickness
    !     if (depth > traverse_depth) then
    !         maxi = i
    !         exit
    !     else
    !         i = i - 1
    !     end if
    ! end do

    do i=N,1,-1
        depth = depth + regolayer(i)%thickness
        depth_diff = depth - traverse_depth
        if(depth_diff > 0) then
            maxi = i
            exit
        end if
    end do


    allocate(poppedarray,source=regolayer(maxi:N))

    !for #1 element of poppedarray, shrink thickness by whatever was lefr over. In corresponding maxi of regolayer, also need to change that.

    poppedarray(1)%thickness = poppedarray(1)%thickness - depth_diff
    regolayer(maxi)%thickness = regolayer(maxi)%thickness - poppedarray(1)%thickness

    ! copy regolayer from 1 to maxi to temp variable, then deallocate regolayer, then movealloc templayer onto regolayer <--may need temp array
    allocate(oldregodata,source=regolayer(1:maxi))
    deallocate(regolayer)
    call move_alloc(oldregodata,regolayer) ! right intents?

    ! if (z <= depth) then
    !     dz = depth - z

    !*****the following lines may still be needed, especially if they deal with thickness:*****

    ! oldregodata                  = regolayer(maxi)
    ! oldregodata%thickness        = z
    ! oldregodata%age(:)           = z / regolayer(maxi)%thickness * regolayer(maxi)%age(:)
    ! recyratio                    = dz / regolayer(maxi)%thickness
    ! regolayer(maxi)%age(:)    = recyratio * regolayer(maxi)%age(:)
    ! regolayer(maxi)%thickness = dz

    !********************


    !call util_push_array(poppedarray,oldregodata) <--not needed; just editing in place
    ! else
    !     z = z - regolayer%thickness
    !     call util_pop_array(regolayer,oldregodata)
    !     call util_push_array(poppedarray,oldregodata)
    !     depth = regolayer%thickness
    ! end if
       
    return
 end subroutine util_traverse_pop_array
 