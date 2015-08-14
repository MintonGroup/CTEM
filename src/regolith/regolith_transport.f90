!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_transport
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Transports material from the inside of the transient crater to
!  the ejecta blanket
!  
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments :
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_transport(user,surf,crater,domain,ejb,ejtble,xp,yp)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_transport
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   type(ejbtype),dimension(:),intent(in)    :: ejb
   integer(I4B),intent(in) :: ejtble
   real(DP),intent(in) :: xp,yp

   ! Internal variables
   real(DP),dimension(4) :: ejecta_corner ! Corners of ejecta block
   real(DP) :: x,y,lrad,loglrad,logtablerad,frac
   integer(I4B) :: i,k
   real(DP),parameter :: maxwellZ = 3.0
   integer(I4B),parameter :: nsteps = 100
   real(DP) :: theta,dtheta

   ! Executable code
   ! Calculate the boundaries of the ejected block that makes this pixel
   do i=1,4
      select case (i)
      case(1)
         x = xp - 0.5 * user%pix
         y = yp - 0.5 * user%pix
      case(2)
         x = xp + 0.5 * user%pix
         y = yp - 0.5 * user%pix
      case(3)
         x = xp - 0.5 * user%pix
         y = yp + 0.5 * user%pix
      case(4)
         x = xp + 0.5 * user%pix
         y = yp + 0.5 * user%pix
      end select

      lrad = sqrt((crater%xl - x)**2 + (crater%yl - yp)**2)
      k = min(int((lrad - crater%frad)/domain%ejbres) + 1,ejtble)
      loglrad=log(lrad)
      logtablerad = log(crater%frad + domain%ejbres*(k-1))
      ! Interpolate back to position inside crater where the corner of the flow came from
      if (k == 1) then
         frac = (ejb(k+1)%erad - ejb(k)%erad)/domain%ejbres
      else
         frac = (ejb(k)%erad - ejb(k-1)%erad)/domain%ejbres
      end if
      ejecta_corner(i) = ejb(k)%erad + frac * (loglrad - logtablerad)
   end do

   ! Now go through the Maxwell-Z streamlines and get the average material
   ! properties in each segment
   theta = 0.5*PI
   dtheta = theta/nsteps
   do i = 1,nsteps

      theta = theta - dtheta
   end do

   return
end subroutine regolith_transport

