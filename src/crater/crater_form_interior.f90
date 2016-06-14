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
subroutine crater_form_interior(user,surfi,crater,lradsq,newelev,melev,&
           thickness_porous_tot,thickness_porous_mare)
   use module_globals
   use module_util
   use module_regolith
   use module_porosity
   use module_crater, EXCEPT_THIS_ONE => crater_form_interior
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: lradsq
   real(DP),intent(in) :: newelev,melev
   real(DP),intent(inout),optional :: thickness_porous_tot, thickness_porous_mare

   ! Internal variables
   real(DP) :: cform,newdem,elchange,pikeD
   integer(I4B) :: layer

   ! Test: internal variables of calculating transient crater's parabola shape
   real(DP) :: trdepth, trparab, trform, trvcorr
   real(DP), parameter :: TRDDRATIO = 0.5_DP !(1.0_DP/3.0_DP + 1.0_DP/4.0_DP) / 2.0_DP 
   real(DP) :: porous_thick
   real(DP) :: x_wall, z_wall, vdiff, cdepth, parabarea, parabside
   type(regodatatype) :: mixedregodata

   ! Executable code

   !change digital elevation map
   cform = crater%vcorr - (crater%parab * lradsq)
   newdem = newelev - cform

   pikeD = 1.044e3_DP * (crater%fcrat * 1e-3_DP)**(0.301_DP) ! Pike (1977)
   if ((crater%fcrat > crater%cxtran * 2) .and. newdem < (melev - pikeD)) then
      newdem = melev - pikeD ! Flatten out the bottom of the crater
   end if
   if (newdem < (melev - user%deplimit)) then
      newdem = melev - user%deplimit ! Flatten out the bottom of the crater
      do layer = 1,user%numlayers ! Remove all pre-existing craters from this current pixel
         call util_remove_from_layer(surfi,layer)
      end do
   end if
   elchange  = newdem - surfi%dem
   surfi%dem = newdem

   if (user%doporosity) then
      call porosity_form_interior(user,surfi,crater,elchange,lradsq,newelev)
   else
      !change ejecta coverage
      surfi%ejcov = max(surfi%ejcov + elchange,0.0_DP)
   end if    
     
   ! Test the parabola shape of a transient crater 
   ! TRDDRATIO: it is similar to DDRATIO, which is the ratio of transient crater's depth to the transient crater diameter
   !            in Jay's cratering book, it is between 1/4 and 1/3. 
   ! trdepth:  it is similar to cdepth, which is the depth of a trasient crater
   !           trdepth = TRDDRATIO * crater%rad * 2.0
   ! trparab:  it is similar to crater%parab, which is the coefficient that determines the shape of a parabola, depending on
   !           the distance of a point on the wall of a crater to the center of a parabola. 
   if (user%doregotrack) then
      call regolith_traverse_pop(elchange,surfi,mixedregodata)
   end if

   !do regotrack: pop stuff out
   !if (user%doregotrack) call regolith_traverse_pop(elchange,surfi)

   return
end subroutine crater_form_interior

