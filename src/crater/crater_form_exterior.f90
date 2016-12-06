!**********************************************************************************************************************************
!
!  Unit Name   : crater_form_exterior
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Forms the crater rim and ejecta blanket on the surface matrix
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
subroutine crater_form_exterior(user,surfi,crater,domain,lradsq,newelev,rd,deltaMi)
   use module_globals
   use module_util
   use module_ejecta
   use module_crater, EXCEPT_THIS_ONE => crater_form_exterior
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(in) :: lradsq,newelev,rd
   real(DP),intent(out) :: deltaMi
   integer(I4B) :: l

   ! Internal variables
   real(DP) :: cform,elchange,belchng,hcorr,lrad

   ! Executable code

   ! interpolate the ejecta blanket thickness
   lrad = sqrt(lradsq)

   ! change digital elevation map
   !cform = (crater%rheight * ((crater%frad / lrad)**RIMDROP))
   cform = (crater%rheight * ((crater%frad / lrad)**rd))

   ! elevation changes
   if (lrad <= crater%frim) then
      hcorr = belchng * ((crater%frim - lrad) / (crater%frim - crater%rad))**20
      elchange = hcorr + cform 
   else 
      elchange = cform 
   end if
   surfi%dem = surfi%dem + elchange
   deltaMi = elchange
   do l = 1, user%numlayers
      surfi%baseline(l) = surfi%baseline(l) + elchange
   end do

   !change ejecta coverage
   !surfi%ejcov = max(surfi%ejcov + elchange,0.0_DP)

   return
end subroutine crater_form_exterior

