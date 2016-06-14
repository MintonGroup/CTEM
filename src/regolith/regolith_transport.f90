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
subroutine regolith_transport(user,surfi,crater,domain,ejb,ejtble,lrad,ebh,comp,popflagi)
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
   INTEGER(I4B),intent(inout) :: popflagi

   ! Internal varialbes
   real(DP) :: melt 
   type(regodatatype) :: newsurfi 

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

   newsurfi%thickness = ebh
   newsurfi%meltfrac = melt
   newsurfi%comp = comp

   call regolith_push(surfi,newsurfi,popflagi)

   return
end subroutine regolith_transport

