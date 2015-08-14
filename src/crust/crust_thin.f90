!**********************************************************************************************************************************
!
!  Unit Name   : crust_thin
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Creates crustal thin area caused by uplifting the mantle  beneath the crater
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : 
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine crust_thin(user,surf,crater,domain,mdepth)
   use module_globals
   use module_crater
   use module_util
   use module_crust, EXCEPT_THIS_ONE => crust_thin
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: mdepth
   
   ! Internal variables
   real(DP)  :: trans_depth ! Transient crater depth
   integer(I4B) :: i,j,xpi,ypi,inc,incsq,iradsq
   real(DP) :: xp,yp,radsq,lradsq,depth

   ! Executable Code
   trans_depth = 2 * THIRD * crater%rad

   ! Don't do anything if the transient crater is not deep enough to cause a  mantle uplift
   if (trans_depth < 0.5_DP * mdepth) return 

   ! determine area to effect
   inc = int(crater%rad/user%pix) !  Maximum distance of crater form
   inc = max(min(inc,PBCLIM*user%gridsize),1)

   radsq = crater%rad**2
   incsq = inc**2

   ! Loop over affected matrix area
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then

            xpi=crater%xlpx+i
            ypi=crater%ylpx+j

            ! Find distance from crater center to current pixel center in real space
            xp = xpi*user%pix
            yp = ypi*user%pix

            lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            depth = surf(xpi,ypi)%dem - surf(xpi,ypi)%mantle

            call crust_cta_generate(user,surf(xpi,ypi),crater,domain,lradsq,depth,mdepth,trans_depth)
         end if
      end do
   end do

   return
   end subroutine crust_thin
