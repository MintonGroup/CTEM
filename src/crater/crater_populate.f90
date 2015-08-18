!**********************************************************************************************************************************
!
!  Unit Name   : crater_populate
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Populates the surface with craters and calls all auxiliary  subroutines
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
subroutine crater_populate(user,surf,crater,domain,prod,vdist,ntrue,vistrue,ntotkilled,truelist,mass,tallycadence,fracdone)
   use module_globals
   use module_seismic
   use module_io
   use module_ejecta
   use module_util
   !use module_crust
   use module_regolith ! simulate regolith reowrking zone
   use module_crater, EXCEPT_THIS_ONE => crater_populate
   implicit none

   ! Arguments
   type(usertype),intent(in)                       :: user
   type(surftype),dimension(:,:),intent(inout)     :: surf
   type(cratertype),intent(inout)                  :: crater
   type(domaintype),intent(inout)                  :: domain
   real(DP),dimension(:,:),intent(in)              :: prod,vdist
   integer(I4B),intent(out)                        :: ntrue
   integer(I4B),intent(out)                        :: vistrue
   integer(I4B),intent(out)                        :: ntotkilled
   real(DP),dimension(:,:),allocatable,intent(out) :: truelist
   real(DP),intent(out)                            :: mass
   integer(I4B),intent(inout)                      :: tallycadence
   real(DP),intent(out)                            :: fracdone

   ! Internal variables
   real(DP)                :: cmin     ! Minimum crater diameter (m)
   real(DP)                :: cmax     ! Maximum crater diameter (m)
   real(DP)                :: imin     ! Minimum impactor diameter (m)
   real(DP)                :: imax     ! Maximum impactor diameter (m)
   real(DP)                :: rmin     ! Depth of smallest crater (m)
   real(DP)                :: rmax     ! Depth of biggest crater (m)
   real(DP)                :: rhmax = -1.0_DP   ! Rim height of biggest crater (m)
   real(DP)                :: rhmin = -1.0_DP    ! Rim height of smallest crater (m)
   real(DP)                :: rhpmax = -1.0_DP   ! Biggest rim height to diameter ratio
   real(DP)                :: rhpmin = -1.0_DP  ! Smallest rim height to diameter ratio
   real(DP)                :: ddmax    ! Maximum depth/diameter 
   real(DP)                :: ddmin    ! Maximum depth/diameter 
   real(DP)                :: melev    ! Mean elevation (m)
   real(DP)                :: xslp,yslp ! Mean slopes
   integer(I8B)            :: icrater  ! Loop counters
   real(DP)                :: mdepth   ! Mean mantle depth below impact site (m)
   integer(I4B)            :: nkilled  ! Number of craters killed in a tally step
   integer(I4B)            :: onum     ! Number of craters observed in a tally step
   integer(I4B)            :: nsincetally ! number of loops since last tally
   integer(I8B)            :: ntotcrat  ! Total number of attempted impacts
   real(DP),dimension(:,:),allocatable  :: tmptruelist
   integer(I4B),parameter  :: TRUECHUNK = 1000000 ! Size of truelist chunks to allocate 
   integer(I4B)            :: truesize
   integer(I4B)            :: craters_since_tally,icrater_last_tally
   integer(I4B)            :: i,j
   real(DP)                :: finterval ! fraction of interval so far completed
   !character(len=MESSAGESIZE) :: message  ! message for the progress bar

   ! ejecta blanket array
   type(ejbtype),dimension(EJBTABSIZE) :: ejb       ! Ejecta blanket lookup table
   integer(I4B) :: ejtble
   ! doregotrack
   real(DP) :: melt
   integer(I4B) :: mratio

   ntotcrat = int(prod(2,domain%smallest_impactor_index),kind=I8B)
   if (user%testflag) then
      ntotcrat = 1
      write(*,*) "Generating a test crater"
      write(*,*) "Dimp = ",user%testimp
      write(*,*) "Vimp = ",user%testvel
      write(*,*) "angle = ",user%testang
      write(*,*) "x offset = ",user%testxoffset
      write(*,*) "y offset = ",user%testyoffset
   else
      write(*,'(" Generating random population of craters. Fewer than: ",I0)') ntotcrat
      write(*,*) "Minimum impactor diameter: ",prod(1,domain%smallest_impactor_index)
      write(*,*) "Minimum crater diameter: ",domain%subcrater_limit
      write(*,*) "Maximum crater diameter: ",domain%biggest_crater
   end if


   ! create crater population
   cmin = domain%biggest_crater
   rmin = 0.0_DP
   cmax = 0.0_DP
   rmax = 0.0_DP
   ntrue = 0
   vistrue = 0
   mass = 0.0_DP
   nsincetally = 0
   ntotkilled = 0
   fracdone = 1.0_DP
   truesize = TRUECHUNK
   allocate(truelist(TRUECOLS,truesize))
   ! begin cratering loop
   if (.not.user%testflag)  then
      pbarival = floor(real(ntrue)/real(PBARRES))
      call io_resetPbar()
   end if
   icrater_last_tally = 0
   icrater = 0
   do while (icrater < ntotcrat)
      icrater = icrater + 1
      pbarpos = ceiling(real(icrater) / real(ntotcrat) * PBARRES)
      ! generate random crater
      call crater_generate(user,crater,domain,prod,vdist,surf)
      if (user%testflag) write(*,*) 'Dcrat = ',crater%fcrat
      if (user%testflag) write(*,*) 'Dtrans = ',crater%rad*2
      if (crater%fcrat > domain%biggest_crater) then ! End the run if the crater is too big
         if (user%killatmaxcrater) then 
            fracdone = real(icrater,kind=DP) / real(ntotcrat,kind=DP)
            write(*,*)
            write(*,'("Ended run at ",F7.2,"% due to crater of size: ",ES13.4)') fracdone * 100,crater%fcrat
            exit
         else
            cycle ! Ignore this big crater
         end if
      end if 

      ! Find the visible crater parameters
      call crater_find_visible(user,crater,domain)

      call ejecta_distance_estimate(user,crater,domain,crater%ejdis) ! Fast but imprecise estimate of the total ejecta distance
                                                                     ! For very steep size distributions, only a fraction of the
                                                                     ! craters are retained. The full ejecta_table_define function
                                                                     ! is very computationally expensive. This function ball-parks
                                                                     ! the total distance to determine if it is worth doing the 
                                                                     ! full calculation later.

      if (((crater%fcrat < domain%smallest_crater) .and. &
          (crater%ejdis < domain%smallest_ejecta)) .or.  &
          (crater%fcrat < domain%subcrater_limit))  cycle ! Ejecta and crater are both too small,so we'll ignore this crater 
      
      ! Crater is big enough to keep, so record it into the true distribution 
      ntrue = ntrue + 1
      if (ntrue > truesize) then  ! Resize the truelist array if necessary
         allocate(tmptruelist(TRUECOLS,truesize))
         tmptruelist = truelist
         deallocate(truelist)
         truesize = truesize + TRUECHUNK
         allocate(truelist(TRUECOLS,truesize))
         truelist(:,1:truesize - TRUECHUNK) = tmptruelist
         deallocate(tmptruelist)
      end if
      truelist(1,ntrue) = crater%fcrat
      truelist(2,ntrue) = crater%imp
      truelist(3,ntrue) = crater%xl
      truelist(4,ntrue) = crater%yl
      truelist(5,ntrue) = crater%impvel
      truelist(6,ntrue) = crater%sinimpang
      mass = mass + crater%impmass

      crater%maxinc = 0
      ! Do seismic shaking
      if (user%doseismic) call seismic_shake(user,surf,crater,domain)

      ! find the average height and slope at crater location
      call crater_averages(user,surf,crater,melev,xslp,yslp,mdepth)
    
      ! Place ejecta onto the surface
      if (crater%ejdis > domain%smallest_ejecta) then ! Estimated size is big enough, so proceed with precise calculation
         if (user%doregotrack) then 
            call ejecta_table_define(user,crater,domain,ejb,ejtble,melt)
            call ejecta_interpolate(crater,domain,crater%frad,ejb(1:ejtble),ejtble,crater%ejrim)
         else 
            call ejecta_table_define(user,crater,domain,ejb,ejtble)
            call ejecta_interpolate(crater,domain,crater%frad,ejb(1:ejtble),ejtble,crater%ejrim)
         end if
         call ejecta_emplace(user,surf,crater,domain,ejb(1:ejtble),ejtble)
      else
         ejtble = 0
      end if

      ! Place crater onto the surface
      if (crater%fcrat > domain%smallest_crater) then
         call crater_emplace(user,surf,crater,domain,melev,xslp,yslp)
         ! Record crater in an available layer as long as it is above the cutoff
         call crater_record(user,surf,crater,melev,xslp,yslp)
         call util_sort_layer(user,surf,crater)
         vistrue = vistrue + 1
         nsincetally = nsincetally + 1
         if (.not.user%testflag) call io_updatePbar("")
      end if


      ! Collapse any remaining unstable slopes
      if (user%docollapse) call crater_slope_collapse(user,surf,crater,domain)

      !if (user%docrustal_thinning) call crust_thin(user,surf,crater,domain,mdepth)


      ! Find out if the current crater is the largest or smallest and if so record it
      if (crater%fcrat > cmax ) then
         imax = crater%imp
         cmax = crater%fcrat
         rhmax = crater%vrim
         if (crater%vcorr <= user%deplimit) then
            rmax = crater%vdepth
         else
            rmax = user%deplimit + crater%vrim 
         end if
         call io_ejecta_table(crater,domain,ejb,ejtble,"ejecta_table_max.dat")
      end if
      if (crater%fcrat < cmin) then
         imin = crater%imp
         cmin = crater%fcrat
         rhmin = crater%vrim
         if (crater%vcorr <= user%deplimit) then
            rmin = crater%vdepth
         else
            rmin = user%deplimit + crater%vrim
         end if
         call io_ejecta_table(crater,domain,ejb,ejtble,"ejecta_table_min.dat")
      end if


      !if (user%doregotrack) then 
      !   mratio = mod( nsincetally, 10 )
      !   if ( mratio == 0 ) then
      !      craters_since_tally = icrater - icrater_last_tally
      !      finterval = craters_since_tally / real(ntotcrat,kind=DP)
      !      call regolith_reworking_zone(user,surf,finterval)
      !   end if
      !end if
   
      if (nsincetally == tallycadence) then
         craters_since_tally = icrater - icrater_last_tally
         finterval = craters_since_tally / real(ntotcrat,kind=DP)
         !if (.not.user%testflag) call ejecta_subcrater_diffusion(user,surf,domain,finterval)
         if (user%doregotrack) then 
         !   call io_write_regodist(user,surf,finterval)
         call regolith_reworking_zone(user,surf,finterval)
         !write(*,*) finterval, craters_since_tally, icrater_last_tally, ntotcrat
         end if
         icrater_last_tally = icrater
         call crater_tally_observed(user,surf,domain,nkilled,onum)
         ntotkilled = ntotkilled + nkilled
         nsincetally = 0
         if (nkilled == 0) then
            tallycadence = tallycadence * 2
         else
            tallycadence = max(1,int(tallycadence * TALLYTARGET / (nkilled / (1.0_DP * user%gridsize**2))))
         end if
      end if

   end do  ! end crater production loop 

   craters_since_tally = icrater - icrater_last_tally
   finterval = craters_since_tally / real(ntotcrat,kind=DP)
   if (ntotcrat == 0) finterval = 1
   if (.not.user%testflag) call ejecta_subcrater_diffusion(user,surf,domain,finterval)

   ! Resize the true crater size array to the actual number of craters produced
   
   ! Display stats
   ddmax = rmax / cmax
   ddmin = rmin / cmin
   rhpmax = rhmax / cmax
   rhpmin = rhmin / cmin
   write(*,*)
   write(*,*) 'Minimum impactor diameter = ',imin
   write(*,*) 'Maximum impactor diameter = ',imax
   write(*,*) 'Minimum crater diameter = ',cmin,' d/D = ',ddmin,' r/D = ', rhpmin
   write(*,*) 'Maximum crater diameter = ',cmax,' d/D = ',ddmax,' r/D = ', rhpmax

   
   return
end subroutine crater_populate

