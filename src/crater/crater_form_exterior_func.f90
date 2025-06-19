!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Forms the exterior raised rim of the crater
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


function crater_form_exterior_func(user,surf,crater,domain,rd,deltaMtot,lastloop) result(ans)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_form_exterior_func
   implicit none

  ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(inout) :: domain
   real(DP),intent(in) :: rd,deltaMtot
   logical,intent(in) :: lastloop
   real(DP) :: ans

   ! Internal variables

   real(DP) :: lradsq,newelev,deltaMp
   integer(I4B) :: xpi,ypi,i,j,k,inc,incsq,iradsq
   real(DP) :: xp,yp,radsq,deltaMi
   type(surftype) :: surfi

   ! Now make the exterior of the crater
   inc = int(crater%frad/user%pix*(domain%small/crater%rimheight)**(-1._DP/RIMDROP)) !  Maximum distance of crater form
   inc = max(min(max(crater%rimdispx,inc),PBCLIM*user%gridsize),1)

   crater%maxinc = max(crater%maxinc,inc)
   incsq = inc**2

   radsq = crater%frad**2

   deltaMp = 0.0_DP

   ! Loop over affected matrix area
   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !$OMP SHARED(inc,radsq,incsq,rd,lastloop) &
   !$OMP SHARED(crater,user,surf) &
   !$OMP REDUCTION(+:deltaMP)
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            ! find elevation and grid point
            newelev = crater%melev + ((i * crater%xslp) + (j * crater%yslp)) * user%pix
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
               

            xp = xpi * user%pix
            yp = ypi * user%pix
            
            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)
            lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

            ! Form interior, rim, and ejecta blanket 
            if (lradsq > radsq) then
               if (lastloop) then
                  call crater_form_exterior(user,surf(xpi,ypi),crater,domain,lradsq,newelev,rd,deltaMi) 
               else
                  surfi = surf(xpi,ypi)
                  call crater_form_exterior(user,surfi,crater,domain,lradsq,newelev,rd,deltaMi) 
               end if
               deltaMp = deltaMp + deltaMi
            end if
         end if
      end do
   end do !end area loopover 
   !$OMP END PARALLEL DO
   ans = deltaMtot + deltaMp
   domain%hmax = maxval(surf(:,:)%dem)
   return

end function crater_form_exterior_func

