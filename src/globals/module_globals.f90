!**********************************************************************************************************************************
!
!  Unit Name   : module_globals
!  Unit Type   : module
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Globally defined constants
!
!  Notes       :  
!
!**********************************************************************************************************************************
module module_globals
implicit none
public

character(len=*),parameter :: CTEMVER = "1.5 DEVELOPMENT"

! Symbolic names for kind types of 4-, 2-, and 1-byte integers:
integer, parameter :: I8B = selected_int_kind(17)
integer, parameter :: I4B = selected_int_kind(9)
integer, parameter :: I2B = selected_int_kind(4)
integer, parameter :: I1B = selected_int_kind(2)

! Symbolic names for kind types of single- and double-precision reals:
integer, parameter :: SP = kind(1.0)
integer, parameter :: DP = kind(1.0d0)

! Frequently used mathematical constants 
real(DP), parameter :: PI      = ACOS(-1.0_DP) 
real(DP), parameter :: SQRT2   = SQRT(2.0_DP) 
real(DP), parameter :: LOGSQRT2 = LOG(SQRT2) 
real(DP), parameter :: SQRT3   = SQRT(3.0_DP) 
real(DP), parameter :: THIRD   = 1.0_DP / 3.0_DP 
real(DP), parameter :: SIXTH   = 1.0_DP / 6.0_DP 
real(DP), parameter :: DEG2RAD = PI / 180.0_DP 


! Maximum string size:
integer(I4B), parameter :: STRMAX = 128

! ASCII character set parameters:
integer(I4B), parameter :: LOWERCASE_BEGIN  = iachar('a')
integer(I4B), parameter :: LOWERCASE_END    = iachar('z')
integer(I4B), parameter :: UPPERCASE_OFFSET = iachar('A') - iachar('a')

! Miscellaneous constants:
real(DP),parameter :: VSMALL  = tiny(1._DP)    ! Very small number
real(DP),parameter :: LOGVSMALL = log(VSMALL)  ! log of a very small number
real(DP),parameter :: VBIG    = huge(1._DP)    ! Very big number
real(DP),parameter :: SMALLFAC = 1e-5_DP       ! Smallest unit of measurement proportional to pixel size
integer(I4B),parameter :: MAXLAYER=20          ! Maximum number of layers (you need roughly 1-2 layers per order of magnitude of 
                                               ! resolution
real(DP),parameter :: TALLYCOVERAGE = 0.01_DP  ! The total area coverage to reach before a tally step is executed
real(DP),parameter :: SUBPIXELCOVERAGE = 0.025_DP ! The total area coverage to reach before a subpixel evaluate step is executed: 0.05_DP
real(DP),parameter :: COOKIESIZE = 3.0_DP      ! Relative size of old crater to new crater that cookie cutting is applied
                                               ! Only craters smaller than COOKIESIZE times the new crater are cookie cut
integer(I2B),parameter :: MAXAGEBINS=60       ! Maximum number of bins in age distribution reset by impact melting
real(DP),parameter :: ALPHA = 0.125_DP
real(DP),parameter  :: DISEJB = 100.0_DP       ! The extent of discontinuous ejecta in the unit of crater radii. It is used in ejecta_table_define.f90
real(DP),parameter :: RCONT = 2.25267_DP       ! Coefficient of continuous ejecta size power law from Moore et al. (1974) - scaled from km to m
real(DP),parameter :: EXPCONT = 1.006_DP       ! Exponent of continuous ejecta size power law from Moore et al. (1974) 
real(DP),parameter :: PF     = 5.0e9          ! The shock pressure exceeding the hungoit elastic limit of common geological materials in Pa. (5 GPa)
real(DP),parameter :: RAD_GP = 20.0_DP         ! The maximum radial position of producing impact glass spherules within a transient crater (unit of crater radii, crater%rad)


type regodatatype 
   real(SP),dimension(MAXAGEBINS) :: age 
   real(DP) :: thickness
   real(DP) :: meltfrac 
   real(DP) :: comp 
   real(DP) :: porosity   ! Porosity: Maximum 1, Minium 0.
   real(DP) :: damage     ! Damage  : Maximum 1, Minium 0.
   real(DP) :: depth      ! Absolute location with respect to the initial surface. 
end type regodatatype
   
type regolisttype
   type(regodatatype) :: regodata
   type(regolisttype),pointer :: next => NULL()
end type

! Derived data type for simulated surface
type surftype
   real(DP),dimension(MAXLAYER) :: diam   ! Crater diameter
   real(SP),dimension(MAXLAYER) :: xl,yl ! Crater center 
   real(SP),dimension(MAXLAYER) :: timestamp ! Crater formation time
   real(DP) :: ejcov                ! Ejecta coverage
   real(DP) :: dem                  ! Digital elevation model
   real(DP) :: mantle               ! Height of mantle (should be smaller than dem)
   type(regolisttype), pointer :: regolayer => null() ! Pointer to the top of the regolith layer stack
   type(regolisttype), pointer :: porolayer => null() ! Pointer to the top of the porosity layer stack
end type surftype

! Derived data type for crater information
type cratertype
   real(DP) :: imp,imprad,impvel,sinimpang,impmass ! Impactor properties
   ! Real domain properties
   real(SP) :: xl,yl               ! Crater center in simulated surface size units 
   real(SP) :: timestamp           ! The formation time of the crater
   real(DP) :: rad                 ! Transient radius
   real(DP) :: grad                ! Strengthless material transient crater radius
   real(DP) :: frad                ! Final crater radius
   real(DP) :: fcrat               ! Final crater diameter
   real(DP) :: ejrad               ! Radius to begin ejecta/raised rim
   real(DP) :: vdepth,vrim         ! parameter for parabolic crater form
   real(DP) :: frim                ! parameter for parabolic crater form
   real(DP) :: rimdis              ! crater form radius (bowl + upturned rim)
   real(DP) :: ejdis               ! ejecta max distance
   real(DP) :: ejrim               ! ejecta height at crater rim
   real(DP) :: cxexp,cxtran        ! simple to complex scaling parameters
   real(DP) :: kdiffterm, saccelterm ! seismic diffusion and accelleration terms
   real(DP) :: continuous          ! Size of the continuous ejecta blanket
   real(DP) :: fe                  ! Equivalent degradation radius
   ! Pixel domain properties
   integer(I4B) :: xlpx,ylpx       ! Crater center in pixels
   integer(I4B) :: fcratpx,fradpx,rimdispx,ejdispx
   integer(I4B) :: maxinc          ! Maximum area affected 
   integer(I4B) :: strflag         ! 0 for regolith, 1 for bedrock
   real(DP) :: melev,xslp,yslp     ! Mean elevation and slopes at pre-existing impact site

   !Crater dimension information  - See Pike (1977) Impact and Explosion Cratering, Fig. 1
   character(STRMAX) :: morphtype ! Type of crater. One of: "SIMPLE", "TRANSITION", "COMPLEX", "PEAKRING", "MULTIRING"
   real(DP) :: rimheight
   real(DP) :: rimwidth 
   real(DP) :: floordepth
   real(DP) :: floordiam
   real(DP) :: peakheight
   integer(I4B),dimension(:),allocatable :: seedarr ! Random number generator seed array
end type cratertype

! Derived data type for domain variables (sizes and dimensions)
type domaintype
   logical      :: initialize ! Is set to true when crater scaling calculations are performed when the domain is initialized
   real(DP)     :: side      ! Length of domain/pseudo-circumference (m)
   real(DP)     :: area      ! Area of domain (m**2)
   real(DP)     :: parea     ! area of a single pixel (m**2)
   real(DP)     :: GM        ! gravitational constant (m**3/s**2)
   real(DP)     :: small     ! Smallest unit of measurement (m)
   integer(I4B) :: smallest_impactor_index  ! Index of the smallest impactor to draw from the production SFD
   integer(I4B) :: smallest_ejecta_index    ! Index of the smallest impactor to produce ejecta otherwise a crater smaller than this will be ignored and treated as a part of vertical mixing model based on the assumption of their insignificant volume of ejecta 
   real(DP)     :: smallest_ejecta_crater
   real(DP)     :: biggest_crater ! Largest crater to generate
   real(DP)     :: smallest_ejecta ! Smallest ejecta to generate (possibly from craters smaller than smallest_produced_crater)
   real(DP)     :: smallest_crater ! Smallest crater that leaves a depression on the surface
   real(DP)     :: subcrater_limit ! Smallest crater that causes any effect on the surface
   real(DP)     :: subpixel_ejecta_thickness ! Average thickness of ejecta produced by subpixel craters
   real(DP)     :: smallest_counted_crater ! Smallest countable crater
   integer(I4B) :: distl     ! Number of bins in the true crater distribution
   integer(I4B) :: pdistl    ! Number of bins in the production distribution
   integer(I4B) :: plo       ! Smallest sqrt(2) index for the lowest bin in the R-plot distribution
   real(DP)     :: ejbres    ! Ejecta blanket lookup table resolution
   integer(I4B) :: pnum      ! size of production function array
   integer(I4B) :: vnum      ! size of velocity distribution array
   integer(I4B) :: rcnum      ! size of real crater list array for quasi-MC
   real(DP)     :: vescsq    ! Escape velocity at target
   integer(I4B) :: vlo       ! Index of lowest valid velocity in the velocity distribution file
   integer(I4B) :: vhi       ! Index of highest valid velocity in the velocity distribution file
   integer(I4B) :: tallycoverage  ! Estimated areal coverage of craters since the last tally
   integer(I4B) :: subpixelcoverage  ! Estimated areal coverage of craters since the last subpixel step
   real(DP)     :: rmsvel            ! Root mean square impact velocity
end type domaintype 

! Derived data type for user input variables
type usertype
   ! Required input variables
   integer(I4B)      :: gridsize  ! Resolution
   integer(I4B)      :: numlayers ! Number of perched layers
   real(DP)          :: pix       ! Pixel size (m)
   real(DP)          :: mu_b      ! Crater scaling exponential constant (ignored for basins)
   real(DP)          :: kv_b      ! Crater scaling linear constant
   real(DP)          :: mu_r      ! Crater scaling exponential constant (ignored for basins)
   real(DP)          :: kv_r      ! Crater scaling linear constant
   integer(I4B)      :: seed      ! Random number generator seed (only used in non-IDL driven mode)
   real(DP)          :: trho_b    ! Target bedrock density
   real(DP)          :: trho_r    ! Target surface regolith layer density
   real(DP)          :: ybar_b    ! Target bedrock strength (Pa)
   real(DP)          :: ybar_r    ! Target regolith strength (Pa)
   real(DP)          :: gaccel    ! Gravitational accel at target
   real(DP)          :: trad      ! Target body radius
   character(STRMAX) :: mat       ! Material type: 1 = silicate, 2 = ice
   real(DP)          :: prho      ! Projectile density
   character(STRMAX) :: sfdfile   ! Name of size distribution file
   character(STRMAX) :: velfile   ! Name of velocity distribution file
   real(DP) :: seisk,cohaccel     ! seismic keff, cohesion breaking acceleration
   
   ! Optional input variables
   logical           :: dorealistic ! Set to T to enable realistic crater morphology. Default is F.
   logical           :: docollapse ! Set T to use the slope collapse model (turning off speeds up the code for testing)
   logical           :: doangle    ! Set to F to only do vertical impacts, otherwise do range of angles (default is T)
   logical           :: doporosity ! Porosity on/off flg. Set to F to turn the model off. Default F. 
   logical           :: doquasimc  ! set to T for quasi-MC run. Default F.
   real(DP)          :: basinimp  ! Impactor size to switch to multiring basin
   real(DP)          :: maxcrat   ! fraction that maximum crater can be relative to grid
   real(DP)          :: deplimit  ! complex crater depth limit

   ! Seismic input variables 
   logical ::  doseismic   ! Set to T if you want to do the seismic shaking model
   real(DP) :: seisq       ! Seismic energy attenuation quality factor (Q)
   real(DP) :: neff        ! impact seismic energy efficiency factor
   real(DP) :: tvel        ! target P-wave (body wave) speed (m/s)
   real(DP) :: tfrac       ! mean free path for seismic wave scattering in medium
   real(DP) :: regcoh      ! target surface regolith layer cohesion

   ! Crater diffusion input parameters
   real(DP) :: Kd1 ! Degradation function coefficient (from Minton et al. (2019))
   real(DP) :: psi ! Degradation function exponent (from Minton et al. (2029))
   real(DP) :: psi2 ! Degradation function large size exponent (from Minton et al. (2020))
   real(DP) :: rbreak ! Degradation function break in exponent (from Minton et al. (2020))
   real(DP) :: fe  ! Scale factor for size of degradation region (from Minton et al. (2019))

   ! Ejecta softening variables
   logical           :: dosoftening  ! Set T to use the extra crater softening model
   real(DP)          :: ejecta_truncation ! Set the number of crater diameters to truncate the ejecta
   logical           :: dorays       ! Set T to use ray model
   logical           :: superdomain  ! Set T to include the superdomain

   ! Regolith tracking variables
   logical           :: doregotrack ! Set T to use the regolith tracking model (EXPERIMENTAL)
   ! Scouring variables
   logical           :: doscour  ! Set T to use the ejecta scouring model (EXPERIMENTAL)
   ! Crustal thinning variables
   logical           :: docrustal_thinning ! Set T to use the crustal thinning model (EXPERIMENTAL)

   logical           :: killatmaxcrater ! Set T to end the run when a crater exceeds the maximum allowable size
   
   ! Test input variables
   logical           :: testflag     ! Set to T if you want a single crater centered on the grid
   real(DP)          :: testimp      ! Test impactor size
   real(DP)          :: testvel      ! Test impactor velocity
   real(DP)          :: testang      ! Test impactor angle
   real(DP)          :: testxoffset  ! Offset of test crater from center in x direction (m)  
   real(DP)          :: testyoffset  ! Offset of test crater from center in y direction (m)   
   logical           :: tallyonly    ! Only run the tally routine (don't generate any new craters)
   logical           :: testtally    ! Set to T to count all non-cookie cut craters, regardless of score
   real(DP)          :: rctime       ! time (in interval units) for emplacement of quasi-MC crater

   ! IDL driver variables
   character(STRMAX) :: impfile      ! Name of impactor size distribution file (impacts per m^2 per y)
   real(DP)          :: interval     ! Length of interval between outputs (y)
   integer(I4B)      :: numintervals ! Total number of intervals
   character(STRMAX) :: runtype      ! Type of run: single or statistical
   logical           :: restart      ! Set to T restart an old run
   logical           :: popupconsole ! Pop up console window every output interval 
   logical           :: saveshaded   ! Output shaded relief images  
   logical           :: saverego     ! Output regolith map images 
   logical           :: savecomp     ! Output composition map images
   logical           :: savepres     ! Output simplified console display images (presentation-compatible images) 
   logical           :: savetruelist ! Save the true cumulative crater distribution for each interval (large file size)
   real(DP)          :: shadedminh   ! Minimum height for shaded relief map (m)
   real(DP)          :: shadedmaxh   ! Maximum height for shaded relief map (m)
   character(STRMAX) :: sfdcompare   ! Type of run: 0 for normal, 1 for statistical (domain is reset between intervals)
   character(STRMAX) :: realcraterlist ! This is only included here so the terminal doesn't return "Unknown parameter"
end type usertype

! Derived data type for the ejecta blanket table elements
type ejbtype
   real(DP) :: lrad     ! Landed radius (m)
   real(DP) :: thick    ! Thickness (m)
   real(DP) :: erad     ! Ejected radius  (m)
   real(DP) :: vesq     ! Ejection velocity squared (m**2 / s**2)
   real(DP) :: angle    ! Ejection angle (deg)
   real(DP) :: meltfrac ! Melt Fraction (melt volume/total ejecta blanket volume at a pixel)
!   real(SP) :: bedrock ! Fraction of bedrock contained in mixture
end type ejbtype

! Progress bar variables
integer(I4B),parameter :: PBARRES = 100
integer(I4B),parameter :: PBARSIZE = 50 
integer(I4B),parameter :: MESSAGESIZE = 32
integer(I4B) :: pbarival
integer(I4B) :: pbarpos
character(len=PBARSIZE) :: pbarchar

! Default file names
character(*),parameter :: DIAMFILE   = 'surface_diam.dat'
character(*),parameter :: TIMEFILE   = 'surface_time.dat'
character(*),parameter :: EJCOVFILE  = 'surface_ejc.dat'
character(*),parameter :: DEMFILE    = 'surface_dem.dat'
character(*),parameter :: REGOFILE   = 'surface_rego.dat'
character(*),parameter :: MELTFILE   = 'surface_melt.dat'
character(*),parameter :: COMPFILE   = 'surface_comp.dat'
character(*),parameter :: STACKNUMFILE = 'surface_stacknum.dat'
character(*),parameter :: AGEFILE = 'surface_age.dat'
character(*),parameter :: STACKPORFILE = 'porosity_stacknum.dat'
character(*),parameter :: POROFILE     = 'porosity_porosity.dat'
character(*),parameter :: DEPTHFILE    = 'porosity_depth.dat'
character(*),parameter :: POROTHICKFILE = 'porosity_thickness.dat'
character(*),parameter :: POSFILE    = 'surface_pos.dat'
character(*),parameter :: TDISTFILE  = 'tdistribution.dat'
character(*),parameter :: TLISTFILE  = 'tcumulative.dat'
character(*),parameter :: ODISTFILE  = 'odistribution.dat'
character(*),parameter :: OLISTFILE  = 'ocumulative.dat'
character(*),parameter :: PDISTFILE  = 'pdistribution.dat'
character(*),parameter :: CRTSCLFILE = 'craterscale.dat'
character(*),parameter :: SFDFILE    = 'production.dat'
character(*),parameter :: USERFILE   = 'ctem.in' 
character(*),parameter :: DATFILE    = 'ctem.dat' !
character(*),parameter :: MASSFILE   = 'impactmass.dat'
character(*),parameter :: RCFILE     = 'craterlist.dat' !not sure if this is where this line should go, but putting it here for now...

! Global variables 
integer(I4B),parameter :: PBCLIM = 1             ! periodic boundary condition limit
integer(I4B),parameter :: SMALLESTCOUNTABLE = 5 ! Minimum pixel diameter for a crater to be considered countable
real(DP),parameter     :: SMALLESTEJECTA = 1.5  ! Minimum number of pixels from center of crater for an ejecta to have any surface effects
integer(I4B),parameter :: TRUECOLS = 7 ! Number of columns in the true crater count array
integer(I4B)           :: NTHREADS = 1 ! Number of OpenMP threads (reset by OpenMP if a parallel environment is detected)
integer(I4B),parameter :: INCPAR = 1   ! Minimum size of inc variables before parallelization kicks in

! Crater scaling parameters
real(DP),parameter :: KT = 0.85_DP             ! Proportionality constant (see Richardson 2009 eqs. 15 & 20)
!real(DP),parameter :: CT = KT * 1.0077158813689795507466256218613060723322903283648264_DP ! KT * (PI*THIRD)**(SIXTH) 
real(DP),parameter :: CT = KT * (PI*THIRD)**(SIXTH) 

real(DP),parameter :: RIMDROP = 4.20_DP        ! Power law index for rim profile 
real(DP), parameter :: EJPROFILE = 3.0_DP      ! Power law index for ejecta profile
real(DP),parameter :: TRSIM = 1.25_DP          ! ?
real(DP),parameter :: EXFAC = 0.1_DP           ! Excavation depth relative to transient crater diameter
real(DP),parameter :: CXEXPS = 1._DP / 0.885_DP - 1.0_DP ! Complex crater scaling exponent (see Croft 1985)
real(DP),parameter :: SIMCOMKS = 16533.8_DP    ! Simple-to-complex transition scaling coefficient for silicate rock
real(DP),parameter :: SIMCOMPS = -1.0303_DP    ! Simple-to-complex transition scaling exponent for silicate rock
real(DP),parameter :: CXEXPI = 0.155_DP        ! ?
real(DP),parameter :: SIMCOMKI = 3081.39_DP    ! Simple-to-complex transition scaling coefficient for ice
real(DP),parameter :: SIMCOMPI = -1.22486_DP   ! Simple-to-complex transition scaling exponent for ice
real(DP),parameter :: SUBPIXFAC = 0.1_DP       ! Subpixel resolution (used for lookup tables and rim creation)
integer(I4B),parameter :: EJBTABSIZE = 1000    ! Lookup table size 
real(DP),parameter :: CRITSLP = 0.7_DP         ! critical slope angle
real(DP),parameter :: COUNTINGRIM = 0.05_DP    ! Fraction inside and outside final diameter to count as rim pixels
real(DP),parameter :: BOWLFRAC = 0.2_DP        ! Fraction of crater interior pixels to use for the bowl-to-rim height calculation
                                               ! (calibrated for Orientale using Potter et al. 2012)
real(DP),parameter :: TRNRATIO = 0.30_DP       ! The ratio of the transient crater depth to the crater diameter.
real(DP),parameter :: KD1PROX = 0.27_DP        ! Intrinsic proximal ejecta degradation function coefficient
real(DP),parameter :: PSIPROX = 2.0_DP         ! Intrinsic proximal ejecta degradation function exponent
real(DP),parameter :: FEPROX  = 1.0_DP         ! Intrinsic proximal ejecta degradation function size scale factor

! Seismic shaking parameters
real(DP),parameter :: SEISFREQ = 20.0_DP    ! seismic wave frequency
real(DP),parameter :: SHEFF = 0.0001_DP    ! linear proportionality constant
real(DP),parameter :: QFAC = 1.80_DP        ! seismic quality factor exponent
real(DP),parameter :: NFAC = 0.294_DP       ! impact efficiency exponent
real(DP),parameter :: PFAC = 0.809_DP       ! impactor diameter exponent
real(DP),parameter :: VFAC = 0.485_DP       ! impactor velocity exponent
real(DP),parameter :: GFAC = 0.487_DP       ! gravitational acceleration exponent
real(DP),parameter :: DFAC = 0.556_DP       ! impact distance exponent

! Crater ray parameters
real(DP),parameter :: rray = 16_DP
integer(I4B),parameter :: Nraymax = 16
real(DP),parameter :: fpeak = 8000_DP ! narrow ray: rw0 propto 1/4
real(DP),parameter :: rayp = 2.0_DP 
integer(I4B),parameter :: rayq = 4
real(DP),parameter :: rayfmult  = (5)**(-4.0_DP / (1.2_DP))  

! quasi-MC test variables
integer(I4B) :: rccount !start line to be read for quasi-MC (should be 0)

end module module_globals
