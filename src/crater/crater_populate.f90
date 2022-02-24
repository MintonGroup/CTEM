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
subroutine crater_populate(user,surf,crater,domain,prod,production_list,vdist,ntrue,vistrue,ntotkilled,truelist,mass, &
                           fracdone,nflux,ntotcrat,curyear,rclist)
   use module_globals
   use module_seismic
   use module_io
   use module_ejecta
   use module_util
   !use module_crust
   use module_regolith ! simulate regolith mixing 
   use module_crater, EXCEPT_THIS_ONE => crater_populate
   implicit none

   ! Arguments
   type(usertype),intent(inout)                    :: user
   type(surftype),dimension(:,:),intent(inout)     :: surf
   type(cratertype),intent(inout)                  :: crater
   type(domaintype),intent(inout)                  :: domain
   real(DP),dimension(:,:),intent(in)              :: prod,vdist
   integer(I8B),dimension(:),intent(inout)         :: production_list
   integer(I4B),intent(out)                        :: ntrue
   integer(I4B),intent(out)                        :: vistrue
   integer(I4B),intent(out)                        :: ntotkilled
   real(DP),dimension(:,:),allocatable,intent(inout) :: truelist
   real(DP),intent(out)                            :: mass
   real(DP),intent(out)                            :: fracdone
   real(DP),dimension(:,:),intent(in)              :: nflux 
   integer(I8B),intent(in)                         :: ntotcrat  ! Total number of attempted impacts
   real(DP),intent(in)                             :: curyear
   real(DP),dimension(:,:), intent(in)             :: rclist !array of 'real' craters for quasiMC

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
   real(DP)                :: timestamp_old
   integer(I8B)            :: icrater  ! Loop counters
   integer(I4B)            :: nkilled  ! Number of craters killed in a tally step
   integer(I4B)            :: onum     ! Number of craters observed in a tally step
   integer(I4B)            :: nsincetally ! number of loops since last tally
   real(DP),dimension(:,:),allocatable  :: tmptruelist
   integer(I4B),parameter  :: TRUECHUNK = 1000000 ! Size of truelist chunks to allocate 
   integer(I4B)            :: truesize
   integer(I4B)            :: craters_since_tally,icrater_last_tally,craters_since_subpixel,icrater_last_subpixel
   integer(I4B)            :: i,j,layer
   real(DP)                :: finterval ! fraction of interval so far completed
   character(len=MESSAGESIZE) :: message  ! message for the progress bar
   real(DP)                :: ejbmass
   logical                 :: makecrater
   real(DP),dimension(user%gridsize,user%gridsize) :: kdiff 
   integer(I4B)            :: oldpbarpos
   real(DP),dimension(:,:),allocatable   :: ejecta_dem
   real(DP)                :: hmax, hmin

   ! ejecta blanket array
   type(ejbtype),dimension(EJBTABSIZE) :: ejb       ! Ejecta blanket lookup table
   integer(I4B) :: ejtble
   ! subpixel utilities
   real(DP),dimension(2,domain%pnum) :: p
   integer(I4B)                      :: craters_since_subpixel_mix, icrater_last_subpixel_mix

   ! doregotrack & age simulation test
   real(DP)              :: melt, clock, age, thick
   real(SP),dimension(user%gridsize, user%gridsize)  :: agetop
   real(SP),dimension(60)                            :: agetot
   type(regolisttype),pointer                        :: current => null()
   real(DP)              :: age_resolution

   if (user%testflag) then
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

   ! read initial quasi-MC position
   if (user%doquasimc) then
      rccount = 1
      user%rctime = rclist(6,rccount)
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
   icrater_last_subpixel = 0
   icrater = 0
   ! Reset age
   clock = 0.0_DP
   finterval = 1.0_DP / real(ntotcrat,kind=DP)
   age     = user%interval
   age_resolution = age / real(MAXAGEBINS)

   ! Reset coverage map
   domain%tallycoverage = 0
   domain%subpixelcoverage = 0
   kdiff = 0.0_DP
   pbarpos = 0
   call io_updatePbar("")
   oldpbarpos = 0
   do while (icrater < ntotcrat)
      makecrater = .true.
      timestamp_old = real(curyear + real(icrater,kind=DP) / real(ntotcrat,kind=DP) * user%interval,kind=SP)
      icrater = icrater + 1
      crater%timestamp = real(curyear + real(icrater,kind=DP) / real(ntotcrat,kind=DP) * user%interval,kind=SP)
      pbarpos = nint(real(icrater) / real(ntotcrat) * PBARRES)
      !if in quasiMC mode: check to see if it's time for a real crater
      if (user%doquasimc) then
         if ((user%rctime > timestamp_old) .and. (user%rctime < crater%timestamp)) then
            write(*,*) "Real crater at ", crater%timestamp
            user%testflag = .true.
            user%testimp = rclist(1, rccount)
            user%testvel = rclist(2, rccount)
            user%testang = rclist(3, rccount)
            user%testxoffset = rclist(4, rccount)
            user%testyoffset = rclist(5, rccount) 
         end if
      end if
      ! generate random crater
      call crater_generate(user,crater,domain,prod,production_list,vdist,surf)
      if (user%testflag) write(*,*) 'Dcrat = ',crater%fcrat
      if (user%testflag) write(*,*) 'Dtrans = ',crater%rad*2
      if (crater%fcrat > domain%biggest_crater) then ! End the run if the crater is too big
         if ( user%testflag .eqv. .false. ) then
            if (user%killatmaxcrater) then 
               fracdone = real(icrater,kind=DP) / real(ntotcrat,kind=DP)
               write(*,*)
               write(*,'("Ended run at ",F7.2,"% due to crater of size: ",ES13.4)') fracdone * 100,crater%fcrat
               exit
            else
               makecrater = .false. ! Ignore this big crater
            end if
         end if
      end if 
      if (user%doquasimc) then
         if (crater%timestamp > user%rctime) then
            user%testflag = .false.
            rccount = rccount + 1
            if (rccount > domain%rcnum) then
               write(*,*) "Real crater list complete."
               user%rctime = 1e30
            else
               user%rctime = rclist(6,rccount)
            end if
         end if
      end if
      if (crater%fcrat > 0.8_DP * domain%side) then
         surf%ejcov  = 0.0_DP
         surf%dem    = 0.0_DP
         do layer = 1,user%numlayers
            surf%diam(layer)   = 0.0_DP
            surf%xl(layer)     = 0.0_SP
            surf%yl(layer)     = 0.0_SP
         end do
         makecrater = .false.
      end if

      if (crater%fcrat < domain%smallest_crater) makecrater = .false.

      if (makecrater) then
         ! Stamp the current time onto the crater
         ! Find the visible crater parameters
         call crater_dimensions(user,crater,domain)

         ! Crater is big enough to keep, so record it into the true distribution 
         ntrue = ntrue + 1
         if (ntrue > truesize) then  ! Resize the truelist array if necessary
            call move_alloc(truelist, tmptruelist)
            truesize = truesize + TRUECHUNK
            allocate(truelist(TRUECOLS,truesize))
            truelist(:,1:truesize - TRUECHUNK) = tmptruelist(:,:)
            deallocate(tmptruelist)
         end if
         truelist(1,ntrue) = crater%fcrat
         truelist(2,ntrue) = crater%imp
         truelist(3,ntrue) = crater%xl
         truelist(4,ntrue) = crater%yl
         truelist(5,ntrue) = crater%impvel
         truelist(6,ntrue) = crater%sinimpang
         truelist(7,ntrue) = crater%timestamp
         mass = mass + crater%impmass

         crater%maxinc = 0
         ! Do seismic shaking
         if (user%doseismic) call seismic_shake(user,surf,crater,domain)
        
         ! find the average height and slope at crater location
         call crater_averages(user,surf,crater)

         ! Place crater onto the surface
         call crater_emplace(user,surf,crater,domain,ejbmass)

         call ejecta_distance_estimate(user,crater,domain,crater%ejdis) ! Fast but imprecise estimate of the total ejecta distance
                                                                        ! For very steep size distributions, only a fraction of the
                                                                        ! craters are retained. The full ejecta_table_define function
                                                                        ! is very computationally expensive. This function ball-parks
                                                                        ! the total distance to determine if it is worth doing the 
                                                                        ! full calculation later.
         ! Place ejecta onto the surface
         if (crater%ejdis > domain%smallest_ejecta) then ! Estimated size is big enough, so proceed with precise calculation
            if (user%doregotrack) then 
               call ejecta_table_define(user,crater,domain,ejb,ejtble,melt)
               !call ejecta_interpolate(crater,domain,crater%frad,ejb(1:ejtble),ejtble,crater%ejrim)
            else 
               call ejecta_table_define(user,crater,domain,ejb,ejtble)
               !call ejecta_interpolate(crater,domain,crater%frad,ejb(1:ejtble),ejtble,crater%ejrim)
            end if
            call ejecta_emplace(user,surf,crater,domain,ejb(1:ejtble),ejtble,ejbmass,age,age_resolution,ejecta_dem)
         else
            ejtble = 0
         end if

         if (user%dorealistic) call crater_realistic_topography(user,surf,crater,domain,ejecta_dem) 
         deallocate(ejecta_dem)

         ! Collapse any remaining unstable slopes
         if (user%docollapse) call crater_slope_collapse(user,surf,crater,domain,(CRITSLP * user%pix)**2,ejbmass)

         ! Record crater in an available layer as long as it is above the cutoff
         call crater_record(user,surf,crater)
         call util_sort_layer(user,surf,crater)
         vistrue = vistrue + 1
         nsincetally = nsincetally + 1
         if (.not.user%testflag) then
            if (pbarpos /= oldpbarpos) then
               call io_updatePbar("")
               oldpbarpos = pbarpos
            end if
         end if

         !if (user%docrustal_thinning) call crust_thin(user,surf,crater,domain,mdepth)
         
         ! Find out if the current crater is the largest or smallest and if so record it
         if (crater%fcrat > cmax ) then
            imax = crater%imp
            cmax = crater%fcrat
            rhmax = crater%vrim
            if (crater%floordepth <= user%deplimit) then
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
            if (crater%floordepth <= user%deplimit) then
               rmin = crater%vdepth
            else
               rmin = user%deplimit + crater%vrim
            end if
            call io_ejecta_table(crater,domain,ejb,ejtble,"ejecta_table_min.dat")
         end if
      end if

      ! Do sub-pixel craters vertical mixing
      ! Do superdomain ray deposits
      finterval = 1.0_DP / real(ntotcrat,kind=DP)
      if (user%doregotrack) then
         call crater_superdomain(user,surf,age,age_resolution,prod,nflux,domain,finterval)
         call regolith_depth_model(user,domain,finterval,nflux,p)
         call regolith_subcrater_mix(user,surf,domain,nflux,finterval,p)
      end if 

      ! Do periodic subpixel processes on the whole grid
      if (.not.user%testflag) then
         if ((domain%subpixelcoverage / real(user%gridsize**2,kind=DP) > SUBPIXELCOVERAGE).or.(icrater == ntotcrat)) then
            domain%subpixelcoverage = 0
            write(message,*) "Subpixel"
            call io_updatePbar(message)
            craters_since_subpixel = icrater - icrater_last_subpixel
            finterval = craters_since_subpixel / real(ntotcrat,kind=DP)
            call crater_subpixel_diffusion(user,surf,nflux,domain,finterval,kdiff)
            icrater_last_subpixel = icrater
         end if
         ! Intermediate tally step 
         if (domain%tallycoverage / real(user%gridsize**2,kind=DP) > TALLYCOVERAGE) then
            domain%tallycoverage = 0
            write(message,*) "Tally"
            call io_updatePbar(message)
            craters_since_tally = icrater - icrater_last_tally
            finterval = craters_since_tally / real(ntotcrat,kind=DP)
            icrater_last_tally = icrater
            call crater_tally_observed(user,surf,domain,nkilled,onum)
            write(message,*) "Tally killed ",nkilled
            call io_updatePbar(message)
            ntotkilled = ntotkilled + nkilled
            nsincetally = 0
         end if
      end if

      hmax = maxval(surf(:,:)%dem)
      hmin = minval(surf(:,:)%dem)
      if (any(surf(:,:)%dem /= surf(:,:)%dem)) then
         write(*,*) crater%imp, crater%impvel, crater%xl, crater%yl, crater%sinimpang
         error stop "Invalid surface elevation detected. Halting."
      end if
   end do  ! end crater production loop 


   call move_alloc(truelist, tmptruelist)
   allocate(truelist(TRUECOLS,ntrue))
   truelist(:,1:ntrue) = tmptruelist(:,1:ntrue)
   deallocate(tmptruelist)
 
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
