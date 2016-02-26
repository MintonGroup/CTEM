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
   real(DP) :: xp,yp,fradsq
   integer(I4B),parameter :: NAVG = 5
   real(DP) :: avgelev
   type(surftype),dimension(NAVG) :: surfavg

   ! Test: porous regime mixing
   real(DP) :: comp_porous, thickness_porous_tot, thickness_porous_mare
   type(regolayertype) :: porouslayer

   ! Executable code

   ! determine area to effect
   inc = int(crater%frad/user%pix*(domain%small/crater%rheight)**(-1._DP/RIMDROP)) !  Maximum distance of crater form
   inc = max(min(max(crater%rimdispx,inc),PBCLIM*user%gridsize),1)

   crater%maxinc = max(crater%maxinc,inc)
   fradsq = crater%frad**2
   incsq = inc**2

   thickness_porous_tot = 0._DP
   thickness_porous_mare = 0._DP

   ! Loop over affected matrix area
   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !$OMP SHARED(inc,fradsq,incsq,melev,xslp,yslp,thickness_porous_tot,thickness_porous_mare) &
   !$OMP SHARED(crater,user,surf) 
   do j=-inc,inc  ! Do the loop in pixel space
      do i=-inc,inc
         ! find distance from crater center
         iradsq = i*i + j*j
         if (iradsq <= incsq) then
            ! find elevation and grid point
            newelev = melev + ((i * xslp) + (j * yslp)) * user%pix
            xpi = crater%xlpx + i
            ypi = crater%ylpx + j
               
            ! periodic boundary conditions
            call util_periodic(xpi,ypi,user%gridsize)

            ! Find distance from crater center to current pixel center in real space
            surfavg = surf(xpi,ypi)
            do k = 1,NAVG
               select case(k)
               case(1)    
                  xp = xpi * user%pix
                  yp = ypi * user%pix
               case(2)
                  xp = (xpi + 0.5_DP) * user%pix 
                  yp = ypi * user%pix
               case(3)
                  xp = (xpi + 0.5_DP) * user%pix 
                  yp = (ypi + 0.5_DP) * user%pix
               case(4)
                  xp = (xpi - 0.5_DP) * user%pix 
                  yp = ypi * user%pix
               case(5)
                  xp = (xpi - 0.5_DP) * user%pix 
                  yp = (ypi - 0.5_DP) * user%pix
               end select 
               lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2

               ! Form interior, rim, and ejecta blanket 
               if (lradsq < fradsq) then 
                  call crater_form_interior(user,surfavg(k),crater,lradsq,newelev,melev,&
                  thickness_porous_tot,thickness_porous_mare)
               else 
                  call crater_form_exterior(user,surfavg(k),crater,domain,lradsq,newelev) 
               end if
            end do
            surf(xpi,ypi) = surfavg(1)
            surf(xpi,ypi)%dem = sum(surfavg%dem) / NAVG
            surf(xpi,ypi)%ejcov = sum(surfavg%ejcov) / NAVG

         end if

      end do
   end do !end area loopover 
   !$OMP END PARALLEL DO
   domain%tallycoverage = domain%tallycoverage + int(fradsq * PI / user%pix**2)

   ! Test: Calculate the theoretical volume difference between transient crater and final crater
   ! dV = PI/4.0 * z_intersect**2 * (R_f/DDRATIO - R_TR/TRDDRATIO)
   ! , where z_intersect is the depth of intersection between transient crater wall
   ! and final crater wall. 
   !comp_porous = thickness_porous_mare / thickness_porous_tot
   !if (user%doregotrack .and. comp_porous == comp_porous .and. thickness_porous_tot > 1.0e-8) then
   !   do j=-inc,inc
   !      do i=-inc,inc
   !         iradsq = i*i + j*j
   !         if (iradsq <= incsq) then
   !            xpi = crater%xlpx + i
   !            ypi = crater%ylpx + j
   !            xp = xpi * user%pix
   !            yp = ypi * user%pix
   !            lradsq = (crater%xl - xp)**2 + (crater%yl - yp)**2
   !            call util_periodic(xpi,ypi,user%gridsize)
   !            if (lradsq < fradsq) then
   !                porouslayer%thickness = surf(xpi,ypi)%regolayer%thickness
   !                porouslayer%comp      = comp_porous
   !                porouslayer%meltfrac  = 0._DP
   !                call regolith_pop(surf(xpi,ypi))
   !                call regolith_push(surf(xpi,ypi),porouslayer)
   !            end if
   !         end if   
   !      end do
   !   end do
   !end if

   return
end subroutine crater_emplace

