!**********************************************************************************************************************************
!
!  Unit Name   : crater_averages
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the average height and slope at the crater location
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
subroutine crater_averages(user,surf,crater,melev,xslp,yslp,mdepth)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_averages
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   type(cratertype),intent(in)  :: crater
   real(DP),intent(out) :: melev,xslp,yslp,mdepth

   ! Internal variables
   integer(I4B) :: mcnt,mx1,mx2,my1,my2
   integer(I4B) :: i,j,inc,incsq,iradsq

   ! Executable code

   mcnt = 0
   melev = 0.0_DP
   xslp = 0.0_DP
   yslp = 0.0_DP
   mdepth = 0.0_DP

   ! determine area to effect
   inc = max(min(crater%fradpx,user%gridsize),1)
   incsq=inc**2

   ! loop over crater area
   do j=-inc,inc 
      do i=-inc,inc

         ! limit to inside of crater radius
         iradsq = i**2 + j**2
         if (iradsq <= incsq) then

            ! define location & distance from crater center
            mx1 = crater%xlpx + i
            my1 = crater%ylpx + j
            mx2 = crater%xlpx + i + 1
            my2 = crater%ylpx + j + 1

            ! periodic boundary conditions
            call util_periodic(mx1,my1,user%gridsize)
            call util_periodic(mx2,my2,user%gridsize)

            ! add pixel element to the total
            melev = melev + surf(mx1,my1)%dem
            mdepth = mdepth + surf(mx1,my1)%dem - surf(mx1,my1)%mantle
            xslp = xslp + ((surf(mx2,my1)%dem - surf(mx1,my1)%dem) / user%pix)
            yslp = yslp + ((surf(mx1,my2)%dem - surf(mx1,my1)%dem) / user%pix)
            mcnt = mcnt + 1
         end if
      end do
   end do

   ! compute working values
   if (mcnt > 1) then 
      melev = melev / real(mcnt,kind=DP)
      xslp = xslp / real(mcnt,kind=DP)
      yslp = yslp / real(mcnt,kind=DP)
      mdepth = mdepth / real(mcnt,kind=DP)
   end if

   return
end subroutine crater_averages

