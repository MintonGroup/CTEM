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
   type(regolayertype),pointer :: current
   type(regolayertype) :: porouslayer
   real(DP) :: z, z0, zmare, ztot
         
   current => surfi%regolayer
   z = current%thickness
   
   if (z < d) then 
      zmare = 0._DP 
      ztot = 0._DP
      
      if (z < d .and. z > VSMALL) then 
          
         do 
          if (.not. associated(current%next)) exit
          if (z <= d) then
             ztot = ztot + current%thickness
             zmare = zmare + current%thickness * current%comp
             current => current%next
             z0 = z
             z = z + current%thickness
          else
             ztot = ztot + (d  - z0)
             zmare = zmare + (d - z0) * current%comp
             exit
          end if
         end do 

         porouslayer%thickness  = ztot
         porouslayer%comp       = zmare / ztot
         porouslayer%meltfrac   = 0.0_DP
         call regolith_traverse_pop(-1.0_DP * ztot, surfi)
         call regolith_push(surfi,porouslayer)
         !write(*,*) '1',surfi%regolayer%thickness, zmare, ztot, zmare/ztot, d
         !if (d > 100.0) write(*,*) ztot, zmare/ztot, d
      end if
   else 
      ztot = d
      zmare = d * surfi%regolayer%comp
      porouslayer%thickness = ztot
      porouslayer%comp      = zmare/ztot
      porouslayer%meltfrac  = 0._DP
      call regolith_traverse_pop(-1.0_DP * ztot, surfi)
      call regolith_push(surfi,porouslayer)
      !write(*,*) '2',surfi%regolayer%thickness, surfi%regolayer%comp, zmare, ztot, zmare/ztot, d
   end if 

   return
end subroutine regolith_mix_porous_regime
