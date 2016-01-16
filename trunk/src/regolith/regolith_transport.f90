!**********************************************************************************************************************************
!
!  Unit Name   : regolith_transport
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Transports material from the inside of the transient crater to
!  the ejecta blanket
!  
!
!  Input
!    Arguments : elchange  ::  elevation change of current ejecta deposited
!                melt      ::  melt fraction of current ejecta deposited 
!                subpixel_ejecta_thickness  ::  the minimum ejecta thickness determined by init_domain.f90  
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_transport(user,surfi,crater,domain,ejb,ejtble,lrad,ebh,comp)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_transport
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(in) :: domain
   integer(I4B),intent(in) :: ejtble
   type(ejbtype),dimension(ejtble),intent(in)   :: ejb
   real(DP),intent(in)          :: lrad,ebh,comp
   !integer(I4B),intent(in)      :: xpi,ypi

   ! Internal varialbes
   real(DP) :: melt 
   type(regolayertype) :: regotop 
   !real(DP) :: minimum_deposition !5.0d-04 average size of impact glass: 500 micronmeter

   ! Melt interpolation variables 
   real(DP)     :: frac,logtablerad,loglrad,logdelta,outeredge,inneredge
   integer(I4B) :: k

   ! Executalbe code

   ! Melt interpolation refered to ejecta_interpolate.f90 

   outeredge = crater%frad + domain%ejbres * (EJBTABSIZE - 0.5_DP)
   inneredge = crater%frad + 0.5_DP * domain%ejbres
   k = max(min(1 + int((lrad - inneredge) / (outeredge - inneredge) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
   loglrad = log(lrad)
   logtablerad = ejb(k)%lrad

   if (k == ejtble) then
      logdelta = logtablerad - ejb(k - 1)%lrad
      frac = (loglrad - ejb(k-1)%lrad) / logdelta
      melt = ejb(k)%meltfrac + ((ejb(k)%meltfrac - ejb(k-1)%meltfrac) * frac)
   else 
      logdelta = ejb(k + 1)%lrad - logtablerad 
      frac = (loglrad - logtablerad) / logdelta 
      melt = ejb(k)%meltfrac - ((ejb(k)%meltfrac - ejb(k+1)%meltfrac) * frac)
   end if 
   regotop%thickness = ebh
   regotop%meltfrac = melt
   regotop%comp = comp

   !minimum_deposition = domain%small * user%gridsize
   !write(*,*) minimum_deposition, domain%small
   !minimum_deposition = 1.0d-03
   !if (ebh >= minimum_deposition) then 
   call regolith_push(surfi,regotop)
   !else if (ebh < minimum_deposition .and. surfi%regolayer%thickness < minimum_deposition) then 
   !        if ( .not. associated(surfi%regolayer%next) ) then 
   !           call regolith_push(surfi,regotop)
   !        else 
   !           surfi%regolayer%meltfrac = ( surfi%regolayer%thickness * surfi%regolayer%meltfrac &
   !                                      + ebh*melt )/(surfi%regolayer%thickness + ebh) 
   !           surfi%regolayer%thickness = surfi%regolayer%thickness + ebh 
   !        end if
   !else if (ebh < minimum_deposition .and. surfi%regolayer%thickness >= minimum_deposition) then
   !        call regolith_push(surfi,regotop)
   !end if

   return
end subroutine regolith_transport

