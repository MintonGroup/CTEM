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

    real(DP) :: maxtdepth, r, tdepth, dz, limit 
    integer(I4B) :: rpix

    !Executable code

    hmeltsheet = vmeltsheet / (nmeltsheet*user%pix*user%pix)
    allocate(newlayer%distvol(1+domain%rcnum))

    rpix = crater%rad / user%pix
    inc = min(rpix,user%gridsize-1)
    limit = 13.0*crater%imprad

    maxtdepth = crater%rad * ((-0.5 * (0.0_DP/crater%rad)**2) + 0.5) ! Parabolic relationship between r and depth of transient crater: 
                                                                    ! depth h = -0.5(r/R)**2 + (1/2) where r is radial distance and R is radius. 
                                                                    ! Max depth of transient crater is assumed to be at r=0, and goes to 0 at r=R.

    do j=-inc,inc
        do i=-inc,inc
            ! find distance from crater center
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
   
            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
   
            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! calculate the depth of the transient crater at this pixel, assuming parabolic shape:
            r = sqrt((crater%xl-xp)**2 + (crater%yl-yp)**2)
            if (r > crater%rad) cycle
            tdepth = -0.5 * (r/crater%rad)**2 + 0.5 ! General use case of parabolic relationship described above

            dz = tdepth * crater%rad ! continuous uplift

            !fill top layer with a breccia lens for megaregolith that is a proportion of crater size (upper limit based on Richardson and Abramov 2020) <--NOTE: Melt sheet has been turned off on this branch and replaced with breccia lens.
            newlayer%ejm = 0.0_DP
            newlayer%thickness = tdepth
            if (newlayer%thickness > 0.684 * crater%floordepth) then
                newlayer%thickness = 0.684 * crater%floordepth
            end if
            newlayer%meltvolume = 0.0_DP !Turn off melt sheet and consider the breccia lens non-melt (mega)regolith
            newlayer%totvolume = newlayer%thickness * (user%pix*user%pix)
            !newlayer%thickness = hmeltsheet
            ! newlayer%meltvolume = vmeltsheet / nmeltsheet
            ! newlayer%totvolume = newlayer%meltvolume
            newlayer%distvol(:) = 0.0_SP
            if(domain%currentqmc) then
                newlayer%distvol(domain%nqmc) = newlayer%meltvolume
            else
                newlayer%distvol(1+domain%rcnum) = newlayer%meltvolume
                newlayer%age(domain%age_counter) = newlayer%meltvolume
            end if
            call util_push_array(surf(xpi,ypi)%regolayer,newlayer)
        end do
    end do


    ! Remove regolith layers in crater wall

    ! do j=-inc,inc
    !     do i=-inc,inc
    !         x_relative = (crater%xl - xp)
    !         y_relative = (crater%yl - yp)
    !         lradsq = x_relative**2 + y_relative**2

    !         if (lradsq > crater%frad**2) cycle
    !         call util_traverse_pop_array(user,surf(xpi,ypi)%regolayer,surf(xpi,ypi)%abselc,poppedarray)
    !         deallocate(poppedarray)

    !         !fill top layer with a breccia lens for megaregolith that is a proportion of crater size (upper limit based on Richardson and Abramov 2020) <--NOTE: Melt sheet has been turned off on this branch and replaced with breccia lens.
    !         newlayer%ejm = 0.0_DP
    !         newlayer%thickness = 0.684 * crater%floordepth !same for all craters since this is an upper limit
    !         newlayer%meltvolume = 0.0_DP !Turn off melt sheet and consider the breccia lens non-melt (mega)regolith
    !         newlayer%totvolume = newlayer%thickness * (user%pix*user%pix)
    !         !newlayer%thickness = hmeltsheet
    !         ! newlayer%meltvolume = vmeltsheet / nmeltsheet
    !         ! newlayer%totvolume = newlayer%meltvolume
    !         newlayer%distvol(:) = 0.0_SP
    !         if(domain%currentqmc) then
    !             newlayer%distvol(domain%nqmc) = newlayer%meltvolume
    !         else
    !             newlayer%distvol(1+domain%rcnum) = newlayer%meltvolume
    !             newlayer%age(domain%age_counter) = newlayer%meltvolume
    !         end if
    !         call util_push_array(surf(xpi,ypi)%regolayer,newlayer)
    !     end do
    ! end do

    deallocate(newlayer%distvol)

    return
end subroutine regolith_interior


