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
   type(domaintype),intent(inout) :: domain
   real(DP),intent(in) :: melev,xslp,yslp

   ! Internal variables
   real(DP) :: lradsq,newelev
   integer(I4B) :: xpi,ypi,i,j,k,inc,incsq,iradsq
   real(DP) :: xp,yp,fradsq,xpii,ypii
   integer(I4B),parameter :: NAVG = 5
   type(surftype),dimension(NAVG) :: surfavg

   ! Executable code

   ! determine area to effect
   inc = int(crater%frad/user%pix*(domain%small/crater%rheight)**(-1._DP/RIMDROP)) !  Maximum distance of crater form
   inc = max(min(max(crater%rimdispx,inc),PBCLIM*user%gridsize),1)

   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   incsq = inc**2

   ! Loop over affected matrix area
   !!$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !!$OMP SHARED(inc,fradsq,incsq,melev,xslp,yslp,thickness_porous_tot,thickness_porous_mare) &
   !!$OMP SHARED(crater,user,surf) 
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
            surfavg = surf(xpi,ypi)

            xpii = crater%xlpx + i
            ypii = crater%ylpx + j
            do k = 1,NAVG ! Average elevation and ejecta coverage over five points inside the grid
               select case(k)
               case(1) 
                  xp = xpii * user%pix
                  yp = ypii * user%pix
               case(2)
                  xp = (xpii + THIRD) * user%pix
                  yp = (ypii + THIRD) * user%pix
               case(3)
                  xp = (xpii + THIRD) * user%pix
                  yp = (ypii - THIRD) * user%pix
               case(4)
                  xp = (xpii - THIRD) * user%pix
                  yp = (ypii + THIRD) * user%pix
               case(5)
                  xp = (xpii - THIRD) * user%pix
                  yp = (ypii - THIRD) * user%pix
               end select
            
               lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

               ! Form interior, rim, and ejecta blanket 
               if (lradsq < fradsq) then 
                  call crater_form_interior(user,surfavg(k),crater,lradsq,newelev,melev)
               else 
                  call crater_form_exterior(user,surfavg(k),crater,domain,lradsq,newelev) 
               end if
            end do
            surf(xpi,ypi) = surfavg(1) ! All properties other than elevation and ejecta coverage are evaluated at the center of the grid
            surf(xpi,ypi)%dem = sum(surfavg%dem) / NAVG
            surf(xpi,ypi)%ejcov = sum(surfavg%ejcov) / NAVG 
         end if
      end do
   end do !end area loopover 
   !!$OMP END PARALLEL DO
   domain%tallycoverage = domain%tallycoverage + int(fradsq * PI / user%pix**2)
   domain%subpixelcoverage = domain%subpixelcoverage + int(fradsq * PI / user%pix**2)

   return
end subroutine crater_emplace

