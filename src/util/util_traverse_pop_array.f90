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
subroutine util_traverse_pop_array(user,regolayer,traverse_depth,poppedarray)
    use module_globals
    use module_util, EXCEPT_THIS_ONE => util_traverse_pop_array
    implicit none
 
    ! Arguments
    !type(regolisttype),pointer   :: regolayer
    type(usertype),intent(in)     :: user
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
    allocate(oldregodata,source=regolayer(maxi:maxi))

    !for #1 element of poppedarray, shrink thickness by whatever was lefr over. In corresponding maxi of regolayer, also need to change that.

    poppedarray(1)%thickness = poppedarray(1)%thickness - depth_diff
    regolayer(maxi)%thickness = regolayer(maxi)%thickness - poppedarray(1)%thickness
    poppedarray(1)%meltvolume = (poppedarray(1)%thickness/oldregodata(1)%thickness) * oldregodata(1)%meltvolume
    poppedarray(1)%distvol(:) = (poppedarray(1)%thickness/oldregodata(1)%thickness) * oldregodata(1)%distvol(:)
    poppedarray(1)%ejm = (poppedarray(1)%thickness/oldregodata(1)%thickness) * oldregodata(1)%ejm
    regolayer(maxi)%meltvolume = (regolayer(maxi)%thickness/oldregodata(1)%thickness) * oldregodata(1)%meltvolume
    regolayer(maxi)%distvol(:) = (regolayer(maxi)%thickness/oldregodata(1)%thickness) * oldregodata(1)%distvol(:)
    regolayer(maxi)%ejm = (regolayer(maxi)%thickness/oldregodata(1)%thickness) * oldregodata(1)%ejm
    poppedarray(1)%age(:) = (poppedarray(1)%thickness/oldregodata(1)%thickness) * oldregodata(1)%age(:)
    regolayer(maxi)%age(:) = (regolayer(maxi)%thickness/oldregodata(1)%thickness) * oldregodata(1)%age(:)

    poppedarray(1)%totvolume = poppedarray(1)%thickness * user%pix * user%pix
    regolayer(maxi)%totvolume = regolayer(maxi)%thickness * user%pix * user%pix

    deallocate(oldregodata)

    ! copy regolayer from 1 to maxi to temp variable, then deallocate regolayer, then movealloc templayer onto regolayer <--may need temp array
    allocate(oldregodata,source=regolayer(1:maxi))
    deallocate(regolayer)
    call move_alloc(oldregodata,regolayer)

       
    return
 end subroutine util_traverse_pop_array
 