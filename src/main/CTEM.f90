!**********************************************************************************************************************************
!
!  Unit Name   : CTEM
!  Unit Type   : program
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : The Cratered Terrain Evolution Model
!
! 
!  Notes       :  This is the main function, which serves mainly to set up the
!  runs and handle input and output. Actually cratering is performed in the 
!  crater_populate subroutine
!
!**********************************************************************************************************************************
program CTEM
use module_globals   
use module_io
use module_init
use module_crater
use module_seismic
use module_ejecta
!$ USE omp_lib
implicit none

! User input variables
type(usertype) :: user

! Domain size variables
type(domaintype) :: domain

! Surface expression grid arrays
type(surftype),dimension(:,:),allocatable :: surf

! Crater properties
type(cratertype) :: crater

! Distribution arrays
real(DP),dimension(:,:),allocatable  :: prod,vdist,pdist,crtscl,truedist,obsdist,truelist
real(DP),dimension(:),allocatable :: obslist
real(SP),dimension(:),allocatable :: original_depth,current_depth,p_score,deviation_sigma
real(SP),dimension(:,:),allocatable :: oposlist

! Miscellaneous variables
character(STRMAX)       :: infile   ! Input file name
logical                 :: restart  ! F = new run (start with a fresh surface)
integer(I8B)            :: totalimpacts ! Total number of impacts ever produced 
integer(I4B)            :: ncount   ! Current count in ctem_driver IDL run
integer(I4B)            :: n        ! Size of random number generator seed array
integer(I4B),dimension(:),allocatable :: seedarr ! Random number generator seed array
real(DP)                :: curyear
real(DP)                :: mass
real(DP)                :: masstot
real(DP)                :: fracdone
integer(I4B)            :: ntrue
integer(I4B)            :: vistrue
integer(I4B)            :: nkilled
integer(I4B)            :: ntotkilled 
integer(I8B)            :: ntotcrat
integer(I4B)            :: onum
!$ real(DP)             :: t1,t2
real(DP),dimension(:,:),allocatable :: nflux

!$ t1 = omp_get_wtime()
call io_splash()
!write(*,*) 'Reading input files'
infile="ctem.in"
call io_input(infile,user)

! Initialize distribution arrays (crater size, number)
!write(*,*) 'Initializing arrays'
call init_dist(user,domain)

allocate(prod(4,domain%pnum))
allocate(crtscl(2,domain%pnum))
allocate(vdist(3,domain%vnum))
allocate(surf(user%gridsize,user%gridsize))

! Read in production impactor population
call io_read_prod(prod,user,domain)

! Read in impactor velocity distribution
call io_read_vdist(vdist,user,domain)

write(*,*) "Initializing simulation domain and determining minimum impactor size"
call init_domain(user,crater,domain,prod,pdist,vdist,crtscl,nflux)

allocate(truedist(6,domain%distl+1))
allocate(obsdist(6,domain%distl+1))

! Reset random number generator
call random_seed
call random_seed(size=n)
allocate(seedarr(n))
call io_read_const(totalimpacts,ncount,curyear,restart,fracdone,masstot,seedarr)
call random_seed(put=seedarr)

! Read in old grid arrays, production function, and velocity distributions
if (restart .or. user%tallyonly) then
   call io_read_surf(user,surf)
else
   call init_surf(user,surf)
end if

if (.not.user%tallyonly) then

   ! Make all the craters!
   if (user%testflag) then
      ntotcrat = 1
   else
      ntotcrat = nint(prod(2,domain%smallest_impactor_index))
   end if
   call crater_populate(user,surf,crater,domain,prod,vdist,ntrue,vistrue,ntotkilled,truelist,mass,fracdone,nflux)

   ! Get the last seed and save it to file
   call random_seed(get=seedarr)
   totalimpacts = totalimpacts + ntotcrat
   call io_write_const(totalimpacts,ncount,curyear,restart,fracdone,masstot,seedarr)
   call crater_tally_true(domain,truelist(:,1:ntrue),ntrue,truedist)
end if

write(*,*) "Tallying craters"
if (.not.user%tallyonly) then
   write(*,*) "Total craters generated:               ",ntotcrat
   write(*,*) "Surface-affecting craters generated:   ",ntrue
   write(*,*) "Visible craters generated:             ",vistrue
end if
call crater_tally_observed(user,surf,domain,nkilled,onum,&
                           obsdist,obslist,oposlist,original_depth,current_depth,deviation_sigma,p_score)
ntotkilled = ntotkilled + nkilled
write(*,*) 'Craters killed during tally: ',ntotkilled
call io_write_tally(truedist,truelist(:,1:ntrue),obsdist,obslist,oposlist,&
                    original_depth,current_depth,deviation_sigma,p_score)
if (.not.user%tallyonly) then
   write(*,*) "Writing surface files"
   call io_write_surf(user,surf)
end if

if (user%testflag) then ! Draw a profile across the crater
   call io_crater_profile(user,surf)
end if
write(*,*) 'Writing output files'

call io_write_dist(pdist,crtscl,domain,mass)

! Deallocate all the allocatables
deallocate(seedarr)
deallocate(surf,prod,vdist,pdist,crtscl,truedist,truelist,obsdist,obslist)
deallocate(oposlist,current_depth,original_depth,p_score,deviation_sigma)

!$ t2 = omp_get_wtime()
!$ write(*,*) 'Timing information'
!$ write(*,*) 'nthreads walltime'
!$ write(*,*) nthreads,t2-t1

end program
