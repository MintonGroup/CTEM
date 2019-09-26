!**********************************************************************************************************************************
!
!  Unit Name   : crater_tally
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Tallies the craters and returns the final distributions
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
subroutine crater_tally_observed(user,surf,domain,nkilled,onum,obsdist,obslist,oposlist,depthdiam)
   use module_globals
   use module_io
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_tally_observed
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(domaintype),intent(in) :: domain
   integer(I4B),intent(out) :: nkilled,onum
   real(DP),dimension(:,:),intent(out),optional  :: obsdist
   real(DP),dimension(:),intent(out),allocatable,optional :: obslist
   real(SP),dimension(:,:),intent(out),allocatable,optional :: oposlist
   real(SP),dimension(:),intent(out),allocatable,optional :: depthdiam

   ! Internal variables
   integer(I4B)                               :: i,j,layer,n,m,craternum,obstot
   type(cratertype) :: crater

   ! master list arrays (these contain information on every pixel that is considered a crater)
   integer(I4B),dimension(user%numlayers*(user%gridsize)**2) :: mlistind
   real(DP),dimension(:),allocatable :: mlist
   real(SP),dimension(:,:),allocatable :: mposlist
   real(SP),dimension(:),allocatable :: tmp_depthdiam
   logical,dimension(:),allocatable :: countable
   real(DP),dimension(:),allocatable :: Kval
   integer(I4B),dimension(:,:),allocatable:: mpxlist
   integer(I4B),dimension(:),allocatable :: ind,mlayerlist
   integer(I4B) :: q,mnum,imnum,maxpix,npix,ntot

   ! tallying arrays (these contain only the number of craters on the grid)
   real(SP),dimension(:,:),allocatable :: poslist
   real(DP),dimension(:),allocatable :: tlist
   real(DP),dimension(:),allocatable :: dis,elev
   integer(I4B),dimension(:),allocatable :: totpix,istart,iend
   integer(I4B) :: tnum ! True number and observed number


   ! Executable code

   ! First build up a master table of every crater for which a pixel still exists on the surface
   mnum=0
   mlistind = 0
   ! Build a mask of all valid crater pixels. This allows us to parallelize the
   ! master list generation
   do n=1,user%gridsize
      do m=1,user%gridsize 
         do layer=1,user%numlayers
            if (surf(m,n)%diam(layer) > 0._DP) then ! Add bottom layer craters to master list
               q = n + user%gridsize * (m - 1) + user%gridsize**2 * (layer - 1)
               mnum = mnum + 1
               mlistind(q) = mnum
            end if
         end do 
      end do
   end do


   allocate(mlist(max(mnum,1)))
   allocate(mposlist(3,max(mnum,1)))
   allocate(mpxlist(2,max(mnum,1)))
   allocate(mlayerlist(max(mnum,1)))
   allocate(ind(max(mnum,1)))
   allocate(istart(max(mnum,1)))
   allocate(iend(max(mnum,1)))

   mlist = 0._DP
   ind = 0

   !$OMP PARALLEL DO DEFAULT(PRIVATE) &
   !$OMP SHARED(user,surf,mlistind) &
   !$OMP SHARED(mlist,mposlist,mpxlist,mlayerlist)
   do n = 1,user%gridsize
      do m = 1,user%gridsize 
         do layer = 1,user%numlayers
            q = n + user%gridsize * (m - 1) + user%gridsize**2 * (layer - 1)
            if (mlistind(q) /= 0) then
               imnum = mlistind(q)
               mlist(imnum) = surf(m,n)%diam(layer) ! The crater diameter, used also as a unique identifier for the crater
               ! Save the crater center coordinates
               mposlist(1,imnum) = surf(m,n)%xl(layer) 
               mposlist(2,imnum) = surf(m,n)%yl(layer) 
               mposlist(3,imnum) = surf(m,n)%timestamp(layer) 
               ! Save the pixel index location
               mpxlist(1,imnum) = m 
               mpxlist(2,imnum) = n 
               mlayerlist(imnum) = layer
            end if
         end do
      end do
   end do 
   !$OMP END PARALLEL DO
   
   ! Sort the master list of all pixels by diameter
   call util_mrgrnk(mlist,ind)

   ! Add up how many unique craters there are
   tnum = 1 
   maxpix = 0
   npix = 1
   istart(1) = 1
   do i=2,mnum
      if (mlist(ind(i)) /= mlist(ind(i-1))) then ! New crater found, skip to next index
         iend(tnum) = i-1
         tnum = tnum + 1
         istart(tnum) = i
         npix = 1
      else
         npix = npix + 1
         if (npix > maxpix) maxpix = npix ! Find out the largest number of pixels occupied by any single crater
      end if
   end do
   iend(tnum) = mnum

   allocate(totpix(tnum))
   allocate(poslist(3,tnum))
   allocate(tlist(tnum))
   allocate(tmp_depthdiam(tnum))
   allocate(countable(tnum))
   allocate(Kval(tnum))
   countable = .false.
   nkilled = 0
   onum = 0
   ! Here we go through the list of craters and take the measures of each one to
   ! be used by the tally subroutine in the next pass
   !$OMP PARALLEL DEFAULT(PRIVATE) IF(tnum > INCPAR) &
   !$OMP SHARED(user,surf) &
   !$OMP SHARED(domain,maxpix,istart,iend,ind,mlist,mposlist,mpxlist,mlayerlist) &
   !$OMP SHARED(tnum,totpix,poslist,tlist,tmp_depthdiam,countable,Kval) &
   !$OMP REDUCTION(+:nkilled) &
   !$OMP REDUCTION(+:onum) 
   !$OMP DO 
   do craternum = 1,tnum
      ! This is the first pixel of this crater, so record the appropriate values
      tlist(craternum) = mlist(ind(istart(craternum)))
      poslist(:,craternum) = mposlist(:,ind(istart(craternum)))
      totpix(craternum) = 1 + iend(craternum) - istart(craternum)
      crater%fcrat = tlist(craternum)
      crater%frad = 0.5_DP * crater%fcrat 
      crater%fradpx = int(crater%frad / user%pix) + 1
      crater%xlpx = int(poslist(1,craternum) / user%pix)
      crater%ylpx = int(poslist(2,craternum) / user%pix)

      call crater_averages(user,surf,crater)

      if ((crater%fcrat < domain%smallest_counted_crater) .or. &
         (totpix(craternum) < int(0.1_DP * 0.25_DP * PI * crater%fcrat / user%pix))) then
         countable(craternum) = .false.
      else 
         Kval(craternum) = crater_get_degradation_state(user,surf,crater,tmp_depthdiam(craternum))
         countable(craternum) = crater_visibility(user,crater,Kval(craternum))
      end if

      if (user%testtally) then ! used for testing the tally system
         countable(craternum) = .true.
      end if

      if (.not.countable(craternum)) then ! Obliterate the crater from the record
         do i=istart(craternum),iend(craternum)
            call util_remove_from_layer(surf(mpxlist(1,ind(i)),mpxlist(2,ind(i))),mlayerlist(ind(i)))
         end do
         nkilled = nkilled + 1
      end if
      if (countable(craternum)) onum = onum + 1
   end do
   !$OMP END DO
   !$OMP END PARALLEL

   ! Bin the observed craters if required 
   if (present(obsdist)) then
      allocate(depthdiam(onum))
      allocate(obslist(onum))
      allocate(oposlist(3,onum))
      ! Reset all the distribution bins
      do i=1,domain%distl
         obsdist(1,i) = 1e3_DP*SQRT2**(domain%plo+(i-1))
         obsdist(2,i) = 1e3_DP*SQRT2**(domain%plo+i)
         obsdist(3,i) = 0.0_DP 
         obsdist(4,i) = 0.0_DP
         obsdist(5,i) = 0.0_DP
         obsdist(6,i) = 0.0_DP
      end do

      ! Count up the observable craters
      obstot = onum
      do craternum=1,tnum
         if (countable(craternum)) then ! Add the crater to the correct differential bin
            ! Find out what bin this crater belongs in
            i = ceiling(log(tlist(craternum)/1e3_DP)/LOGSQRT2) - domain%plo 
            obsdist(3,i) = obsdist(3,i) + log(tlist(craternum)) ! Geometric mean diameter (intermediate step)
            obsdist(4,i) = obsdist(4,i) + 1 ! Differential number
            obslist(obstot) = tlist(craternum)
            oposlist(:,obstot) = poslist(:,craternum)
            depthdiam(obstot) = tmp_depthdiam(craternum)
            obstot = obstot - 1
         end if
      end do

      ! Finalize the bin values
      do i = 1,domain%distl
         if (obsdist(4,i) > 0._DP) then
            obsdist(3,i) = exp(obsdist(3,i) / obsdist(4,i)) ! Geometric mean diameter (final step)
         else
            obsdist(3,i) = sqrt(obsdist(1,i) * obsdist(2,i))
         end if
         obsdist(5,i) = sum(obsdist(4,i:domain%distl)) ! Cumulative number
         obsdist(6,i) = (obsdist(4,i)*obsdist(3,i)**3)/(domain%area*(obsdist(2,i)-obsdist(1,i))) ! R-value
      end do
   end if

   deallocate(poslist)
   deallocate(tlist)
   deallocate(totpix)
   deallocate(istart)
   deallocate(iend)
   deallocate(mlist)
   deallocate(mposlist)
   deallocate(mpxlist)
   deallocate(mlayerlist)
   deallocate(tmp_depthdiam)
   deallocate(ind)
   deallocate(countable)
   deallocate(Kval)

   return
end subroutine crater_tally_observed

