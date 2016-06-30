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
subroutine crater_emplace(user,surf,crater,domain,melev,xslp,yslp,ejbmass)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_emplace
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(inout) :: domain
   real(DP),intent(in) :: melev,xslp,yslp,ejbmass

   ! Internal variables
   real(DP) :: lradsq,newelev
   integer(I4B) :: xpi,ypi,i,j,k,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq,xpii,ypii
   real(DP) :: deltaMtot,deltaMi,deltaMp,rd
   type(surftype) :: surfi
   real(DP),parameter :: KMAX = 100 ! Number of iteratios of the external calculation to converge on a mass-conserving rim profile
   real(DP),parameter :: TOL = 1.0e-10_DP ! Tolerance on the ratio of displaced volume prior to and after crater_form_exterior call
   logical :: lastloop 

   ! Executable code

   ! determine area to effect
   ! First make the interior of the crater
   inc = max(min(crater%rimdispx,PBCLIM*user%gridsize),1) + 1
   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   deltaMtot = ejbmass
   ! This loop may not be paralelizable because of the linked list operation inside crater_form_interior
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find distance from crater center
         
         ! find elevation and grid point
         newelev = melev + ((i * xslp) + (j * yslp)) * user%pix
         xpi = crater%xlpx + i
         ypi = crater%ylpx + j
            

         xp = xpi * user%pix
         yp = ypi * user%pix
         
         ! periodic boundary conditions
         call util_periodic(xpi,ypi,user%gridsize)

         lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2
         if (lradsq <= fradsq) then
            call crater_form_interior(user,surf(xpi,ypi),crater,lradsq,newelev,melev,deltaMi)
            deltaMtot = deltaMtot + deltaMi
         end if
      end do
   end do


   ! Now make the exterior of the crater
   inc = int(crater%frad/user%pix*(domain%small/crater%rheight)**(-1._DP/RIMDROP)) !  Maximum distance of crater form
   inc = max(min(max(crater%rimdispx,inc),PBCLIM*user%gridsize),1)

   crater%maxinc = max(crater%maxinc,inc)
   incsq = inc**2

   rd = RIMDROP ! initial guess of the raised rim profile
   
   if (abs(deltaMtot) > 0.0_DP) then ! Avoid a divide by zero condition in the iteration
      lastloop = .false.
   else
      lastloop = .true.
   end if

   ! Iterate until we converge on a volume-conserving solution to the rim profile
   do k = 1,KMAX
      deltaMp = 0.0_DP
      ! Loop over affected matrix area
      !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
      !$OMP SHARED(inc,fradsq,incsq,melev,xslp,yslp,k,rd,lastloop) &
      !$OMP SHARED(crater,user,surf) &
      !$OMP REDUCTION(+:deltaMP)
      do j=-inc,inc  ! Do the loop in pixel space
         do i=-inc,inc
            ! find distance from crater center
            iradsq = i*i + j*j
            if (iradsq <= incsq) then
               ! find elevation and grid point
               newelev = melev + ((i * xslp) + (j * yslp)) * user%pix
               xpi = crater%xlpx + i
               ypi = crater%ylpx + j
                  

               xp = xpi * user%pix
               yp = ypi * user%pix
               
               ! periodic boundary conditions
               call util_periodic(xpi,ypi,user%gridsize)
               lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

               ! Form interior, rim, and ejecta blanket 
               if (lradsq > fradsq) then
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
      if (lastloop) exit
      if (abs(1.0_DP - abs(deltaMp / deltaMtot)) < TOL) lastloop = .true.
      rd = 2.0_DP - (2.0_DP - rd) * abs(deltaMp / deltaMtot) ! adust the rim drop exponent until volume convergence is reached
      if (rd /= rd) then 
         rd = RIMDROP
         lastloop = .true.
      end if
   end do
   domain%tallycoverage = domain%tallycoverage + int(fradsq * PI / user%pix**2)
   domain%subpixelcoverage = domain%subpixelcoverage + int(fradsq * PI / user%pix**2)

   return
end subroutine crater_emplace

