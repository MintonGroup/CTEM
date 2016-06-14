!**********************************************************************************************************************************
!
!  Unit Name   : regolith_mix_porous_regime
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Simulate mixing driven by cratering mechanics inside a porous regime which is between 
!                the floor of a final crater and transient crater. 
!
!  Input
!    Arguments : 
!
!  Output
!    Arguments : surf : Surface expression matrix
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_mix_porous_regime(user,surfi,d)
   use module_globals
   use module_util
   use module_regolith, EXCEPT_THIS_ONE => regolith_mix_porous_regime
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   real(DP),intent(in) :: d

   ! Regotrack Internals
   type(regodatartype) :: porouslayer
   real(DP) :: z, z0, ztot, zmare
         
   z = surfi%regolayer%thickness
   ztot  = 0._DP
   zmare = 0._DP
   z0    = 0._DP

   do 
    if (.not. associated(surfi%regolayer%next)) exit
    if (z <= d) then
       ztot  = ztot  + surfi%regolayer%thickness
       zmare = zmare + surfi%regolayer%thickness * surfi%regolayer%comp
       call util_pop(surfi)
       z0 = z
       z = z + surfi%regolayer%thickness
    else
       ztot  = ztot  + (d - z0)
       zmare = zmare + (d - z0) * surfi%regolayer%comp
       call regolith_traverse_pop(-1.0_DP * (d - z0), surfi)
       exit
    end if
   end do 

   porouslayer%thickness  = ztot
   porouslayer%comp       = zmare / ztot
   porouslayer%meltfrac   = 0.0_DP
   call util_push(surfi, porouslayer)

   return
end subroutine regolith_mix_porous_regime
