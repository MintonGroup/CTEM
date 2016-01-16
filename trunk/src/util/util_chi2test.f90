!**********************************************************************************************************************************
!
!  Unit Name   : util_chi2
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculates the reduced chi**2 between two arrays
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
function util_chi2test(data1,data2,nbins,p) result(ans)
use module_globals
use module_util, EXCEPT_THIS_ONE => util_chi2test
implicit none
interface gamq
   function gammq(a,x) result(ans)
   use module_globals
   implicit none
   real(DP),intent(in) :: a,x
   real(DP) :: ans
   end function gammq
end interface

integer(I4B),intent(in) :: nbins
real(DP),dimension(nbins),intent(in) :: data1,data2
real(DP),intent(out) :: p
real(DP) :: ans
real(DP) :: DoF
integer(I4B) :: ite

ans=0.0_DP
DoF=nbins
do ite=1,nbins
   if ((data1(ite)==0.0_DP).and.(data2(ite)==0.0_DP)) then
      DoF=DoF-1._DP
   else
      ans=ans+(data1(ite)-data2(ite))**2/(data1(ite)+data2(ite))
   end if
end do
p=gammq(0.5_DP*DoF,0.5_DP*ans)
return
end function util_chi2test

function gammq(a,x) result(ans)
use module_globals
implicit none
interface gser
   subroutine gser(gamser,a,x,gln)
   use module_globals
   implicit none
   real(DP),intent(out) :: gamser,gln
   real(DP),intent(in) :: a,x
   end subroutine gser
end interface gser
interface gcf
   subroutine gcf(gammcf,a,x,gln)
   use module_globals
   implicit none
   real(DP),intent(out) :: gammcf,gln
   real(DP),intent(in) :: a,x
   end subroutine gcf
end interface gcf
real(DP),intent(in) :: a,x
real(DP) :: ans
real(DP) :: gammcf,gamser,gln

if (x<a+1._DP) then ! Series
   call gser(gamser,a,x,gln)
   ans=1._DP-gamser
else
   call gcf(gammcf,a,x,gln)
   ans=gammcf
end if
return
end function gammq

subroutine gser(gamser,a,x,gln)
use module_globals
implicit none
real(DP),intent(out) :: gamser,gln
real(DP),intent(in) :: a,x
integer(I4B),parameter :: itmax=100
real(DP),parameter :: eps=3e-7_DP
integer(I4B) :: n
real(DP) :: ap,del,summ

gln=log(gamma(a))
if (x<=0.) then
   gamser=0._DP
   return
end if
ap=a
summ=1._DP/a
del=summ
do n=1,itmax
   ap=ap+1
   del=del*x/ap
   summ=summ+del
   if (abs(del)<abs(summ)*eps) exit
end do
gamser=summ*exp(-x+a*log(x)-gln)
return
end subroutine gser

subroutine gcf(gammcf,a,x,gln)
use module_globals
implicit none
real(DP),intent(out) :: gammcf,gln
real(DP),intent(in) :: a,x
integer(I4B),parameter :: itmax=100
real(DP),parameter :: eps=3e-7_DP
real(DP),parameter :: fpmin=1e-30_DP
integer(I4B) :: i
real(DP) :: an,b,c,d,del,h

gln=log(gamma(a))
b=x+1._DP-a
c=1._DP/fpmin
d=1._DP/b
h=d
do i=1,itmax
   an=-i*(i-a)
   b=b+2._DP
   d=an*d+b
   if (abs(d)<fpmin) d=fpmin
   c=b+an/c
   if (abs(c)<fpmin) c=fpmin
   d=1._DP/d
   del=d*c
   h=h*del
   if (abs(del-1._DP)<eps) exit
end do
gammcf=exp(-x+a*log(x)-gln)*h
return
end subroutine gcf

