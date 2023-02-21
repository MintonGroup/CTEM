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
!    Arguments : user,surf,domain,crater,inc,nmeltsheet,vmeltsheet
!
!  Output
!    Arguments : surf
!           
! 
!  Notes       : nmeltsheet is the number of pixels that will contain the melt sheet; vmeltsheet is the melt sheet volume
!
!**********************************************************************************************************************************
subroutine regolith_interior(user,surf,crater,domain,incval,nmeltsheet,vmeltsheet)
    use module_globals 
    use module_util
    use module_regolith, EXCEPT_THIS_ONE => regolith_interior
    implicit none

    !Arguments
    type(usertype),intent(in) :: user
    type(surftype),dimension(:,:),intent(inout) :: surf
    type(domaintype),intent(in) :: domain
    type(cratertype),intent(in) :: crater
    integer(I4B),intent(in)     :: incval, nmeltsheet
    real(DP),intent(in)         :: vmeltsheet

    !internal variables
    integer(I4B) xpi,ypi,i,j,inc,incsq,iradsq
    real(DP) :: lradsq, x_relative, y_relative, xp, yp, hmeltsheet 
    type(regodatatype),dimension(:),allocatable :: poppedarray
    type(regodatatype) :: newlayer

    !Executable code

    hmeltsheet = vmeltsheet / (nmeltsheet*user%gridsize*user%gridsize)
    allocate(newlayer%meltdist(domain%rcnum))

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

            !fill top layer with melt sheet of given thickness hmeltsheet
            newlayer%meltfrac = 1.0_DP
            newlayer%ejm = 0.0_DP
            newlayer%ejmf = 0.0_DP
            newlayer%thickness = hmeltsheet
            newlayer%meltvolume = vmeltsheet / nmeltsheet
            newlayer%totvolume = newlayer%meltvolume
            newlayer%meltdist(:) = 0.0_SP
            if(domain%currentqmc) then
                newlayer%meltdist(domain%nqmc) = newlayer%meltfrac
            end if
            call util_push_array(surf(xpi,ypi)%regolayer,newlayer)
        end do
    end do

    deallocate(newlayer%meltdist)

    return
end subroutine regolith_interior


