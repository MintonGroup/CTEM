!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_interpolate
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Interpolate ejecta blanket thickness 
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
subroutine ejecta_interpolate(crater,domain,lrad,ejb,ejtble,ebh,vsq,theta,erad,melt)
   use module_globals
   use module_util
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_interpolate
   implicit none

   ! Arguments
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in)  :: lrad
   integer(I4B),intent(in) :: ejtble
   type(ejbtype),dimension(ejtble),intent(in) :: ejb
   real(DP),intent(out) :: ebh
   real(DP),intent(out),optional :: vsq,theta,erad
   real(DP),intent(out),optional :: melt

   ! Internals
   real(DP)     :: frac,logtablerad,loglrad,logdelta,outeredge,inneredge
   integer(I4B) :: k

   ! Executable code

   ! Locate ourselves in the table
   outeredge = crater%frad + domain%ejbres * EJBTABSIZE
   inneredge = crater%frad 
   k = max(min(1 + int((lrad - inneredge) / (outeredge - inneredge) * (EJBTABSIZE - 1.0_DP)),ejtble),1)
   loglrad = log(lrad)
   logtablerad = ejb(k)%lrad 

   ! Interpolate in logspace (this saves on the number of table elements we need)
   if (k == ejtble) then
      logdelta = logtablerad - ejb(k - 1)%lrad
      frac = (loglrad - logtablerad) / logdelta
      ebh = ejb(k)%thick - ((ejb(k)%thick - LOGVSMALL) * frac)
      if (present(vsq)) vsq = ejb(k)%vesq - ((ejb(k)%vesq - LOGVSMALL) * frac)
      if (present(theta)) theta = ejb(k)%angle - ((ejb(k)%angle - LOGVSMALL) * frac)
      if (present(melt)) melt = ejb(k)%meltfrac - ((ejb(k)%meltfrac - LOGVSMALL) * frac)
      if (present(erad)) erad = ejb(k)%erad - ((ejb(k)%erad - LOGVSMALL) * frac)
   else
      logdelta = ejb(k + 1)%lrad - logtablerad 
      frac = (loglrad - logtablerad) / logdelta
      ebh = ejb(k)%thick - ((ejb(k)%thick - ejb(k + 1)%thick) * frac)
      if (present(vsq)) vsq = ejb(k)%vesq - ((ejb(k)%vesq - ejb(k + 1)%vesq) * frac)
      if (present(theta)) theta= ejb(k)%angle - ((ejb(k)%angle - ejb(k + 1)%angle) * frac)
      if (present(melt)) melt = ejb(k)%meltfrac - ((ejb(k)%meltfrac - ejb(k+1)%meltfrac) * frac) 
      if (present(erad)) erad = ejb(k)%erad - ((ejb(k)%erad - ejb(k+1)%erad) * frac) 
   end if
   ebh = exp(ebh) 
   if (lrad > crater%ejdis) ebh = 0._DP
  
   return

end subroutine ejecta_interpolate
