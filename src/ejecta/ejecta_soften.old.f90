!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_soften
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Softens the terrain under the ejecta  using a box filter model where the size of the box is proportional to the 
!                thickness of the ejecta  
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
subroutine ejecta_soften(user,surf,inc,indarray,cumulative_elchange)
   use module_globals
   use module_util
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_soften
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   integer(I4B),intent(in) :: inc
   integer(I4B),dimension(2,-inc:inc,-inc:inc),intent(in) :: indarray
   real(DP),dimension(-inc:inc,-inc:inc),intent(inout) :: cumulative_elchange 

   ! Internal variables
   integer(I4B) :: xpi,ypi,i,j,boxpx,boxi,boxj,deltai,deltaj,boxpxsq,boxrempx,boxrempxsq,boxradsq,bnum,remnum,ifilt
   real(DP),parameter :: filter_factor = 10.0_DP
   real(DP) :: demavg,boxsize,boxrem,demremavg,fac,ebhmax,boxmax,ebh
   logical :: wrapcheck

   !$OMP PARALLEL DO DEFAULT(PRIVATE) &
   !$OMP SHARED(inc,cumulative_elchange,indarray,ebh) &
   !$OMP SHARED(user,surf) 
   do j=-inc,inc
      do i=-inc,inc
         ebh = cumulative_elchange(i,j)
         boxsize = (filter_factor * ebh)/user%pix
         boxrempx = ceiling(boxsize)
         boxpx = boxrempx - 1 
         boxrempxsq = boxrempx**2
         boxpxsq = boxpx**2
         boxrem = 1.0_DP - (boxrempx - boxsize)
         
         demavg = 0.0_DP 
         if (boxrempx <= 0) then
            cycle
         end if
         demremavg = 0.0_DP
         ! Draw a box around the current pixel of a size determined by the
         ! filter factor
         bnum = 0
         remnum = 0
         !if ((i+boxrempx) > inc .or. (i-boxrempx < -inc) .or. (j+boxrempx) > inc .or. (j-boxrempx) < -inc) then
            wrapcheck = .true.
         !else
         !   wrapcheck = .false.
         !end if
         do boxj=-boxrempx,boxrempx
            do boxi=-boxrempx,boxrempx
               boxradsq = boxi*boxi + boxj*boxj
               if (boxradsq <= boxrempxsq) then
                  if (wrapcheck) then
                     xpi = indarray(1,i,j) + boxi
                     ypi = indarray(2,i,j) + boxj
                     call util_periodic(xpi,ypi,user%gridsize)
                  else  
                     xpi = indarray(1,i+boxi,j+boxj)
                     ypi = indarray(2,i+boxi,j+boxj)
                  end if
                  if ((xpi < 1) .or. (xpi > user%gridsize) .or. (ypi < 1) .or.  (ypi > user%gridsize)) then
                     write(*,*) "xpi,ypi out of bounds!",xpi,ypi
                  end if

                  demremavg = demremavg + surf(xpi,ypi)%dem
                  remnum = remnum + 1
                  if (boxradsq <= boxpxsq) then
                     demavg = demavg + surf(xpi,ypi)%dem
                     bnum = bnum + 1
                  end if
               end if
            end do
         end do
         demavg = demavg / bnum
         demremavg = demremavg  / remnum

         ! This "feathers" the edges of the box to prevent artificial stepping at integer values of the filter size
         demavg = (demremavg - demavg) * boxrem + demavg

         ! Replace the existing DEM value with the average DEM at this location
         xpi = indarray(1,i,j)
         ypi = indarray(2,i,j)
         cumulative_elchange(i,j) = cumulative_elchange(i,j) - surf(xpi,ypi)%dem + demavg
      end do
   end do
   !$OMP END PARALLEL DO

return
end subroutine ejecta_soften

