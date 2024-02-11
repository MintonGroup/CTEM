!**********************************************************************************************************************************
!
!  Unit Name   : util_sort_layer
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Sorts the craters down in the layer
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
subroutine util_sort_layer(user,surf,crater)
   use module_globals
   use module_util, EXCEPT_THIS_ONE => util_sort_layer
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(in) :: crater

   ! Internal variables
   integer(I4B),dimension(user%numlayers) :: ind
   real(DP),dimension(user%numlayers) :: tempdiam
   real(SP),dimension(user%numlayers) :: tempxpos,tempypos,temptime
   integer(I4B) :: i,j,k,inc,incsq,mx,my,iradsq

   inc = min(crater%rimdispx,(user%gridsize - 1)/2)
   incsq = inc*inc


   ! Executable code
   !$OMP PARALLEL DO DEFAULT(PRIVATE) IF(inc > INCPAR) &
   !$OMP SHARED(user,crater,surf) &
   !$OMP SHARED(inc,incsq)
   do j=-inc,inc
      do i=-inc,inc
         iradsq = i*i+j*j
         if (iradsq <= incsq) then
            mx = crater%xlpx + i
            my = crater%ylpx + j
            call util_periodic(mx,my,user%gridsize)
           
            ! Temporarily store layer data
            tempdiam = surf(mx,my)%diam(1:user%numlayers)
            tempxpos = surf(mx,my)%xl(1:user%numlayers)
            tempypos = surf(mx,my)%yl(1:user%numlayers)
            temptime = surf(mx,my)%timestamp(1:user%numlayers)

            ! Sort the layers by crater diameter
            call util_mrgrnk(tempdiam,ind)

            do k=1,user%numlayers
               surf(mx,my)%diam(k)=tempdiam(ind(k))
               surf(mx,my)%xl(k)=tempxpos(ind(k))
               surf(mx,my)%yl(k)=tempypos(ind(k))
               surf(mx,my)%timestamp(k)=temptime(ind(k))
            end do
         end if
      end do
   end do
   !$OMP END PARALLEL DO

   return
end subroutine util_sort_layer

