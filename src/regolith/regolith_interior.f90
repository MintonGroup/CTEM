!**********************************************************************************************************************************
!
!  Unit Name   : regolith_interior
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Pops off and destroys layers in the crater inteorior and fills with melt sheet
!  
!
!  Input
!    Arguments : user,surf,crater,inc
!
!  Output
!    Arguments : surf
!           
! 
!  Notes       :
!
!**********************************************************************************************************************************
subroutine regolith_interior(user,surf,crater,incval)
    use module_globals 
    use module_util
    use module_regolith, EXCEPT_THIS_ONE => regolith_interior
    implicit none

    !Arguments
    type(usertype),intent(in) :: user
    type(surftype),dimension(:,:),intent(inout) :: surf
    type(cratertype),intent(in) :: crater
    integer(I4B),intent(in)     :: incval

    !internal variables
    integer(I4B) xpi,ypi,i,j,inc,incsq,iradsq
    real(DP) :: lradsq, x_relative, y_relative, xp, yp 
    type(regodatatype),dimension(:),allocatable :: poppedarray

    !Executable code

    inc = incval

    do j=-inc,inc
        do i=-inc,inc
            iradsq = i**2 + j**2

            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            xp = xpi * user%pix
            yp = ypi * user%pix

            call util_periodic(xpi,ypi,user%gridsize)
            x_relative = (crater%xl - xp)
            y_relative = (crater%yl - yp)
            lradsq = x_relative**2 + y_relative**2

            if (lradsq > crater%frad**2) cycle
            call util_traverse_pop_array(surf(xpi,ypi)%regolayer,surf(xpi,ypi)%abselc,poppedarray)
            deallocate(poppedarray)

            !Adding the melt sheet goes here! :)
        end do
    end do

    return
end subroutine regolith_interior


