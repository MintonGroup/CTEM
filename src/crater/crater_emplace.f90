!**********************************************************************************************************************************
!
!  Unit Name   : crater_emplace
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds the visible crater parabolic parameters, rim, and  rim upturn distance
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
subroutine crater_emplace(user,surf,crater,domain,melev,xslp,yslp)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_emplace
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: melev,xslp,yslp

   ! Internal variables
   real(DP) :: lradsq,newelev
   integer(I4B) :: xpi,ypi,i,j,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq

   ! Executable code

   ! determine area to effect
   inc = int(crater%frad/user%pix*(domain%small/crater%rheight)**(-1._DP/RIMDROP)) !  Maximum distance of crater form
   inc = max(min(max(crater%rimdispx,inc),PBCLIM*user%gridsize),1)

   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   incsq = inc**2

   ! Loop over affected matrix area
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            ! find elevation and grid point
            newelev = melev + ((i * xslp) + (j * yslp)) * user%pix
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j

            ! Find distance from crater center to current pixel center in real space
            xp = xpi * user%pix
            yp = ypi * user%pix
            
            lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! Form interior, rim, and ejecta blanket 
            if (lradsq < fradsq) then 
               call crater_form_interior(user,surf(xpi,ypi),crater,lradsq,newelev,melev)
            else 
               call crater_form_exterior(user,surf(xpi,ypi),crater,domain,lradsq,newelev) 
            end if

         end if

      end do
   end do !end area loopover 

   return
end subroutine crater_emplace

