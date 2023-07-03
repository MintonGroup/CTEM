!**********************************************************************************************************************************
!
!  Unit Name   : crater_form_interior
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Finds the visible crater parabolic parameters, rim, and  rim upturn distance
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
subroutine crater_form_interior(user,surfi,crater,x_relative, y_relative ,newelev,deltaMi)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_form_interior
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(inout) :: crater
   real(DP),intent(in) :: x_relative, y_relative 
   real(DP),intent(in) :: newelev
   real(DP),intent(out) :: deltaMi

   ! Internal variables
   real(DP) :: cform,newdem,elchange,r
   integer(I4B) :: layer

   ! An array for popped data 
   !type(regolisttype),pointer :: poppedlist
   type(regodatatype),dimension(:),allocatable :: poppedarray


   ! Empirical crater shape parameters from Fassett et al. (2014)
   real(DP),parameter :: r_floor = 0.2_DP
   real(DP),parameter :: simple_depth_diam = 0.181_DP
   real(DP),parameter :: r_rim = 0.98_DP
   real(DP),parameter :: inner_c0 = -0.229_DP 
   real(DP),parameter :: inner_c1 =  0.228_DP 
   real(DP),parameter :: inner_c2 =  0.083_DP 
   real(DP),parameter :: inner_c3 = -0.039_DP

   real(DP),parameter :: outer_c0 =  0.188_DP
   real(DP),parameter :: outer_c1 = -0.187_DP
   real(DP),parameter :: outer_c2 =  0.018_DP 
   real(DP),parameter :: outer_c3 =  0.015_DP

   !real(DP) :: c0,c1,c2,c3,flrad,rh,fld


   ! Executable code
   r = sqrt(x_relative**2+y_relative**2) / crater%frad

   cform = crater_profile(user,crater,r)
   newdem = newelev + cform 

   if (newdem < (crater%melev - crater%floordepth)) then 
      newdem = crater%melev - crater%floordepth ! Flatten out the bottom of the crater regardless of the local slope
      do layer = 1,user%numlayers ! Remove all pre-existing craters from the flat floor
         call util_remove_from_layer(surfi,layer)
      end do
   end if

   newdem = min(newdem,surfi%dem) ! Only allow excavation, no deposition
   elchange  = newdem - surfi%dem
   deltaMi = elchange
   surfi%dem = newdem
   surfi%abselc = abs(elchange)
   
   !change ejecta coverage
   surfi%ejcov = max(surfi%ejcov + elchange,0.0_DP)

   return
end subroutine crater_form_interior

