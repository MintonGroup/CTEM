!**********************************************************************************************************************************
!
!  Unit Name   : io_input
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Read in input file
!
!  Input
!    Arguments : infile : input filename
!
!  Output
!    Arguments : 
! 
!  Notes       :  Subroutine sets several global variables based on user input file
!
!**********************************************************************************************************************************

subroutine io_input(infile,user)
   use module_globals
   use module_util
   use module_io, EXCEPT_THIS_ONE => io_input
   implicit none

   ! Arguments
   character(*), intent(in)  :: infile
   type(usertype),intent(inout) :: user

   ! Internals
   integer(I4B), parameter :: LUN = 7
   integer(I4B)            :: ierr, ilength, ifirst, ilast,i
   character(STRMAX)       :: line, token
   integer(I4B), parameter :: numrequired=18
   character(STRMAX),dimension(numrequired),parameter :: requiredvar = (/"GRIDSIZE ",&
                                                          "NUMLAYERS", &
                                                          "PIX      ", &
                                                          "SEED     ", &
                                                          "GACCEL   ", &
                                                          "TRAD     ", &
                                                          "MU_R     ", &
                                                          "KV_R     ", &
                                                          "YBAR_R   ", &
                                                          "TRHO_R   ", &
                                                          "MU_B     ", &
                                                          "KV_B     ", &
                                                          "YBAR_B   ", &
                                                          "TRHO_B   ", &
                                                          "MAT      ", &
                                                          "PRHO     ", &
                                                          "SFDFILE  ", &
                                                          "VELFILE  "/)

   integer(I4B), parameter :: seismic_numrequired=5
   character(STRMAX),dimension(seismic_numrequired),parameter :: seismic_requiredvar = (/ "SEISQ  ", &
                                                                          "NEFF   ", &
                                                                          "TVEL   ", &
                                                                          "TFRAC  ", &
                                                                          "REGCOH " /)

   logical,dimension(numrequired) :: ismissing=.true.
   logical,dimension(seismic_numrequired) :: seismic_ismissing=.true.
   logical :: haltflag
   ! Executable code
   
   ! Set up default values for optional variables
   user%deplimit = huge(0._DP)
   user%testflag = .false.
   user%testimp = 50._DP
   user%testvel = 17e3_DP
   user%testang = 90._DP
   user%doseismic = .false.
   user%docollapse = .true.
   user%testxoffset = 0._DP
   user%testyoffset = 0._DP
   user%docrustal_thinning = .false.
   user%doscour = .false.
   user%tallyonly = .false.
   user%dosoftening = .true.
   user%diffusion_const = 0.0_DP
   user%doregotrack = .false.
   user%basinimp = huge(0._DP)
   user%maxcrat = 1.00_DP
   user%doangle = .true.
   user%testtally = .false.
   user%killatmaxcrater = .false.
   user%countingmodel = "FASSETT"
   user%tallystart = 5000
   
   open(unit=LUN,file=infile,status="old",iostat=ierr)
   if (ierr /= 0) then
   write(*,*) "Unable to open file ",trim(infile)
      stop
   end if
100 format(A)
   do
      read(LUN, 100, iostat = ierr, end = 1) line
      line = adjustl(line)
      ilength = len_trim(line)
      if ((ilength /= 0) .and. (line(1:1) /= "!")) then
         ifirst = 1
         call io_get_token(line, ilength, ifirst, ilast, ierr)
         token = line(ifirst:ilast)
         call util_toupper(token)
         select case (token)
         case (trim(requiredvar(1)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%gridsize
            ismissing(1)=.false.
         case (trim(requiredvar(2)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%numlayers
            ismissing(2)=.false.
         case (trim(requiredvar(3)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%pix   
            ismissing(3)=.false.
         case (trim(requiredvar(4)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%seed
            ismissing(4)=.false.
         case (trim(requiredvar(5)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%gaccel 
            ismissing(5)=.false.
         case (trim(requiredvar(6)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%trad 
            ismissing(6)=.false.
         case (trim(requiredvar(7)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%mu_r
            ismissing(7)=.false.
         case (trim(requiredvar(8)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%kv_r
            ismissing(8)=.false.
         case (trim(requiredvar(9)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%ybar_r
            ismissing(9)=.false.
         case (trim(requiredvar(10)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%trho_r
            ismissing(10)=.false.
         case (trim(requiredvar(11)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%mu_b
            ismissing(11)=.false.
         case (trim(requiredvar(12)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%kv_b
            ismissing(12)=.false.
         case (trim(requiredvar(13)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%ybar_b 
            ismissing(13)=.false.
         case (trim(requiredvar(14)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%trho_b
            ismissing(14)=.false.
         case (trim(requiredvar(15)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%mat
            call util_toupper(user%mat) 
            ismissing(15)=.false.
         case (trim(requiredvar(16)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%prho 
            ismissing(16)=.false.
         case (trim(requiredvar(17)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%sfdfile
            ismissing(17)=.false.
         case (trim(requiredvar(18)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token,*) user%velfile
            ismissing(18)=.false.
         case ("COUNTINGMODEL")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%countingmodel
            call util_toupper(user%countingmodel)
         case ("DEPLIMIT")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%deplimit
         case ("TESTFLAG")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testflag
         case ("TESTIMP")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testimp
         case ("TESTVEL")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testvel
         case ("TESTANG")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testang
         case ("DOSEISMIC")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%doseismic
         case ("DOCOLLAPSE")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%docollapse
         case ("TESTXOFFSET")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testxoffset
         case ("TESTYOFFSET")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testyoffset
         case ("DOCRUSTAL_THINNING")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%docrustal_thinning
         case ("DOSCOUR")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%doscour
         case ("DOSOFTENING")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%dosoftening
         case ("DIFFUSION_CONST")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%diffusion_const
         case ("DOREGOTRACK")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%doregotrack
         case ("TALLYONLY")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%tallyonly
         case ("BASINIMP")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%basinimp
         case ("MAXCRAT")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%maxcrat
         case ("DOANGLE")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%doangle
         case ("TESTTALLY")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%testtally
         ! Porosity model
         case ("POROSITYFLG")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%porosityflg
         ! Seismic variables (only required if doseismic is set to .true.
         case ("KILLATMAXCRATER")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%killatmaxcrater
         ! Seismic variables (only required if doseismic is set to .true.
         case (trim(seismic_requiredvar(1)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%seisq 
            seismic_ismissing(1)=.false.
         case (trim(seismic_requiredvar(2)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%neff 
            seismic_ismissing(2)=.false.
         case (trim(seismic_requiredvar(3)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%tvel
            seismic_ismissing(3)=.false.
         case (trim(seismic_requiredvar(4)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%tfrac
            seismic_ismissing(4)=.false.
         case (trim(seismic_requiredvar(5)))
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%regcoh
            seismic_ismissing(5)=.false.
         ! IDL inputs
         case ("IMPFILE")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%impfile
         case ("INTERVAL")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%interval
         case ("NUMINTERVALS")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%numintervals
         case ("RUNTYPE")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%runtype
            call util_toupper(user%runtype)
         case ("RESTART")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%restart
         case ("POPUPCONSOLE")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%popupconsole
         case ("SAVESHADED")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%saveshaded
         case ("SAVEREGO")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%saverego
         case ("SAVECOMP")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%savecomp
         case ("SAVEPRES")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%savepres
         case ("SAVETRUELIST")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%savetruelist
         case ("SFDCOMPARE")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%sfdcompare
         case ("SHADEDMINH")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%shadedminh
         case ("SHADEDMAXH")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%shadedmaxh
         case ("TALLYSTART")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%tallystart
         !**************************************************************************
         ! The following is for backwards compatibility with older style input files
         !**************************************************************************
         case("DENS")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%trho_b
            ismissing(14)=.false.
            write(*,*) "Warning! Outdated input file detected!"
            write(*,*) "Please change input variable DENS to TRHO_B."
         case ("REGRHO")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%trho_r
            ismissing(10)=.false.
            write(*,*) "Warning! Outdated input file detected!"
            write(*,*) "Please change input variable REGRHO to TRHO_R."
         case ("YBAR")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%ybar_b 
            ismissing(13)=.false.
            write(*,*) "Warning! Outdated input file detected!"
            write(*,*) "Please change input variable YBAR to YBAR_B."
         case ("RYBAR")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%ybar_r
            ismissing(9)=.false.
            write(*,*) "Warning! Outdated input file detected!"
            write(*,*) "Please change input variable RYBAR to YBAR_R."
         case ("MU")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%mu_b
            user%mu_r = user%mu_b
            ismissing(7)=.false.
            ismissing(11)=.false.
            write(*,*) "Warning! Outdated input file detected!"
            write(*,*) "Please change input variable MU to MU_R and MU_B."
         case ("KV")
            ifirst = ilast + 1
            call io_get_token(line, ilength, ifirst, ilast, ierr)
            token = line(ifirst:ilast)
            read(token, *) user%kv_b
            user%kv_r = user%kv_b
            ismissing(8)=.false.
            ismissing(12)=.false.
            write(*,*) "Warning! Outdated input file detected!"
            write(*,*) "Please change input variable KV to KV_R and KV_B."
         !**************************************************************************
         case default
            write(*, 100, advance = "no") "Unknown parameter -> "
            write(*, *) token
         end select
      end if
   end do
1  close(LUN)
   haltflag=.false.
   do i=1,numrequired
      if (ismissing(i)) then 
         write(*,*) 'Required variable ',trim(adjustl(requiredvar(i))),' is missing.'
         haltflag=.true.
      end if
   end do

   if (user%doseismic) then
      do i=1,seismic_numrequired
         if (seismic_ismissing(i)) then 
            write(*,*) 'Required variable ',trim(adjustl(seismic_requiredvar(i))),' is missing.'
            haltflag=.true.
         end if
      end do
   end if

   !**************************************************************************
   ! The following is for backwards compatibility with older style input files
   !**************************************************************************
   select case(user%mat)
   case("1")
      write(*,*) "Warning! Outdated input file detected!"
      write(*,*) "Please change input variable MAT from 1 to ROCK"
      user%mat = "ROCK"
   case("2")
      write(*,*) "Warning! Outdated input file detected!"
      write(*,*) "Please change input variable MAT from 2 to ICE"
      user%mat = "ICE"
   end select
   !**************************************************************************

   select case(user%countingmodel)
   case("FASSETT")
   case default
      write(*,*) 'Unknown counting model ',trim(adjustl(user%countingmodel)),' specified. Using FASSETT instead.'
      user%countingmodel = "FASSETT"
   end select

   if (user%numlayers > MAXLAYER) then
      write(*,*) 'NUMLAYERS is too high. Must not exceed ',MAXLAYER
      stop
   end if

   if (haltflag) then
      write(*,*) 'Execution stopped because of missing required variables.'
      stop
   end if

   return

end subroutine io_input
