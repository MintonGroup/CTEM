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
subroutine crater_tally_observed(user,surf,domain,nkilled,onum,obsdist,obslist,oposlist,&
                        original_depth,current_depth,deviation_sigma,p)
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
   real(SP),dimension(:),intent(out),allocatable,optional :: original_depth,current_depth,p,deviation_sigma

   ! Internal variables
   integer(I4B)                               :: i,j,layer,n,m,craternum,obstot

   ! master list arrays (these contain information on every pixel that is considered a crater)
   integer(I4B),dimension(user%numlayers*(user%gridsize)**2) :: mlistind
   real(DP),dimension(:),allocatable :: mlist,melevation
   integer(I2B),dimension(:),allocatable :: misrim
   real(SP),dimension(:,:),allocatable :: mposlist
   real(SP),dimension(:),allocatable :: moriginal_depth_list,mbaseline
   real(SP),dimension(:),allocatable :: tmp_original_depth,tmp_current_depth
   real(SP),dimension(:),allocatable :: tmp_p,tmp_deviation_sigma
   logical,dimension(:),allocatable :: countable
   integer(I4B),dimension(:,:),allocatable:: mpxlist
   integer(I4B),dimension(:),allocatable :: ind,mlayerlist
   integer(I4B) :: q,mnum,imnum,maxpix,npix,ntot

   ! tallying arrays (these contain only the number of craters on the grid)
   real(SP),dimension(:,:),allocatable :: poslist
   real(DP),dimension(:),allocatable :: tlist
   real(DP),dimension(:),allocatable :: dis,elev
   integer(I4B),dimension(:),allocatable :: elevind,totpix,istart,iend
   integer(I4B) :: tnum ! True number and observed number
   integer(I4B) :: inc,xpi,ypi
   real(DP) :: avgrim,avgbowldepth,avgbowldepth_orig,xk,Mkm1,xkdd,Mkm1dd
   real(DP) :: xkdev,Mkm1dev,deviation,sigma
   logical :: killable
   integer(I4B) :: nrim,nbowl
   integer(I4B) :: nprof,Ni
   real(DP) :: totavg,avgelev,profres,avgdis,rim,bowl,rad

   ! Executable code


   profres = 0.5_DP * user%pix

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
   allocate(misrim(max(mnum,1)))
   allocate(melevation(max(mnum,1)))
   allocate(mbaseline(max(mnum,1)))
   allocate(mposlist(2,max(mnum,1)))
   allocate(mpxlist(2,max(mnum,1)))
   allocate(moriginal_depth_list(max(mnum,1)))
   allocate(mlayerlist(max(mnum,1)))
   allocate(ind(max(mnum,1)))
   allocate(istart(max(mnum,1)))
   allocate(iend(max(mnum,1)))

   mlist = 0._DP
   ind = 0
   misrim = 0

   !$OMP PARALLEL DO DEFAULT(PRIVATE) &
   !$OMP SHARED(user,surf,mlistind) &
   !$OMP SHARED(mlist,misrim,melevation,mbaseline,mposlist,mpxlist,moriginal_depth_list,mlayerlist)
   do n = 1,user%gridsize
      do m = 1,user%gridsize 
         do layer = 1,user%numlayers
            q = n + user%gridsize * (m - 1) + user%gridsize**2 * (layer - 1)
            if (mlistind(q) /= 0) then
               imnum = mlistind(q)
               mlist(imnum) = surf(m,n)%diam(layer) ! The crater diameter, used also as a unique identifier for the crater
               misrim(imnum) = surf(m,n)%isrim(layer) ! Is this pixel part of the original rim? 1 = yes, 0 = no
               melevation(imnum) = surf(m,n)%dem ! Save the elevation of this pixel
               ! Save the crater center coordinates
               mposlist(1,imnum) = surf(m,n)%xl(layer) 
               mposlist(2,imnum) = surf(m,n)%yl(layer) 
               ! Save the pixel index location
               mpxlist(1,imnum) = m 
               mpxlist(2,imnum) = n 
               moriginal_depth_list(imnum) = surf(m,n)%original_depth(layer) ! The original elevation of the crater
               mbaseline(imnum) = surf(m,n)%baseline(layer)
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
   allocate(poslist(2,tnum))
   allocate(tlist(tnum))
   allocate(tmp_original_depth(tnum))
   allocate(tmp_current_depth(tnum))
   allocate(tmp_p(tnum))
   allocate(tmp_deviation_sigma(tnum))
   allocate(countable(tnum))
   countable = .false.
   nkilled = 0
   onum = 0
   !$OMP PARALLEL DEFAULT(PRIVATE) IF(tnum > INCPAR) &
   !$OMP SHARED(user,surf) &
   !$OMP SHARED(maxpix,istart,iend,ind,mlist,mposlist,misrim,melevation,mbaseline,moriginal_depth_list,mpxlist,mlayerlist) &
   !$OMP SHARED(tnum,totpix,poslist,tlist,tmp_original_depth,tmp_current_depth,tmp_p,tmp_deviation_sigma,countable) &
   !$OMP REDUCTION(+:nkilled) &
   !$OMP REDUCTION(+:onum) 
   !allocate(dis(4*user%gridsize**2))
   !allocate(elev(4*user%gridsize**2))
   !allocate(elevind(4*user%gridsize**2))
   ! Here we go through the list of craters and take the measures of each one to
   ! be used by the tally subroutine in the next pass
   !$OMP DO 
   do craternum=1,tnum
      ! This is the first pixel of this crater, so record the appropriate values
      tlist(craternum) = mlist(ind(istart(craternum)))
      poslist(:,craternum) = mposlist(:,ind(istart(craternum)))
      totpix(craternum) = 1 + iend(craternum) - istart(craternum)
      if (totpix(craternum) < int(0.1_DP * 0.25_DP * PI * tlist(craternum) / user%pix)) then
         countable(craternum) = .false.
         killable = .true.
         tmp_p(craternum) = 0.0_DP
      else
         nrim = 0
         nbowl = 0
         nprof = 0
         bowl = 0.0_DP
         rim = 0.0_DP
         inc = ceiling(0.5_DP * tlist(craternum) / user%pix * 1.25_DP)  
         do j = -inc, inc
            do i = -inc, inc
               rad = sqrt((i**2 + j**2)*1._DP) * user%pix / tlist(craternum)
               if (rad < 0.1_DP) then
                  xpi = int(poslist(1,craternum) / user%pix) + i
                  ypi = int(poslist(2,craternum) / user%pix) + j
                  call util_periodic(xpi,ypi,user%gridsize)
                  bowl = bowl + surf(xpi,ypi)%dem
                  nbowl = nbowl + 1 
               else if (rad > 0.4_DP) then
                  xpi = int(poslist(1,craternum) / user%pix) + i
                  ypi = int(poslist(2,craternum) / user%pix) + j
                  call util_periodic(xpi,ypi,user%gridsize)
                  rim = rim + surf(xpi,ypi)%dem
                  nrim = nrim + 1 
               end if
            end do
         end do
         rim = rim / nrim 
         bowl = bowl / nbowl 
         tmp_p(craternum) = (rim - bowl) / tlist(craternum)

!               if ((i**2 + j**2) < inc**2) then
!                  xpi = int(poslist(1,craternum) / user%pix) + i
!                  ypi = int(poslist(2,craternum) / user%pix) + j
!                  call util_periodic(xpi,ypi,user%gridsize)
!                  nprof = nprof + 1
!                  elev(nprof) = surf(xpi,ypi)%dem
!                  dis(nprof) = sqrt((i**2 + j**2)*1.d0) * user%pix
!                  totavg = totavg + elev(nprof)
!               end if
!            end do
!         end do
!         totavg = totavg / nprof
!         elev = elev - totavg
!         call util_mrgrnk(dis(1:nprof),elevind(1:nprof))
!         avgelev = 0.0_DP
!         avgdis = 0.0_DP
!         Ni = 1
!         ntot = 0
!         i = 0
!         rim = 0.0_DP
!         bowl = 0.0_DP
!         nrim = 0
!         nbowl = 0
!         outer: do 
!            do 
!               i = i + 1
!               if (i > nprof) exit outer
!               if (dis(elevind(i)) > profres * Ni) then
!                  avgdis = avgdis / ntot / tlist(craternum)
!                  avgelev = avgelev / ntot / tlist(craternum)
!                  if (avgdis >= 0.4_DP) then
!                     rim = rim + avgelev
!                     nrim = nrim+ 1
!                  else if (avgdis <= 0.1_DP) then
!                     bowl = bowl + avgelev
!                     nbowl = nbowl + 1
!                  end if
!                  avgelev = elev(elevind(i))
!                  avgdis = dis(elevind(i))
!                  ntot = 1
!                  exit
!               else
!                  avgelev = avgelev + elev(elevind(i))
!                  avgdis = avgdis + dis(elevind(i))
!                  ntot = ntot + 1
!               end if
!            end do
!            Ni = Ni + 1
!         end do  outer


         if (tmp_p(craternum) > 0.05_DP) then
            countable(craternum) = .true.
            killable = .false.
         else  
            countable(craternum) = .false.
            killable = .true.
         end if   
      end if
      ! TESTING
      if (user%testtally) then
         countable = .true.
         killable = .false.
      end if
      if (killable) then ! Obliterate the crater from the record
         do i=istart(craternum),iend(craternum)
            call util_remove_from_layer(surf(mpxlist(1,ind(i)),mpxlist(2,ind(i))),mlayerlist(ind(i)))
         end do
         nkilled = nkilled + 1
      end if
      if (countable(craternum)) onum = onum + 1
      !deallocate(csurf)
   end do
   !$OMP END DO
   !deallocate(dis)
   !deallocate(elev)
   !$OMP END PARALLEL

   ! Bin the observed craters if required 
   if (present(obsdist)) then
      allocate(original_depth(onum))
      allocate(current_depth(onum))
      allocate(p(onum))
      allocate(deviation_sigma(onum))
      allocate(obslist(onum))
      allocate(oposlist(2,onum))
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
            original_depth(obstot) = tmp_original_depth(craternum)
            current_depth(obstot) = tmp_current_depth(craternum)
            p(obstot) = tmp_p(craternum)
            deviation_sigma(obstot) = tmp_deviation_sigma(craternum)
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
   deallocate(misrim)
   deallocate(melevation)
   deallocate(mbaseline)
   deallocate(mposlist)
   deallocate(mpxlist)
   deallocate(moriginal_depth_list)
   deallocate(mlayerlist)
   deallocate(tmp_original_depth)
   deallocate(tmp_current_depth)
   deallocate(tmp_p)
   deallocate(tmp_deviation_sigma)
   deallocate(ind)
   deallocate(countable)

   return
end subroutine crater_tally_observed

