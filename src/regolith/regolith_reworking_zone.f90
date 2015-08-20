!**********************************************************************************************************************************
!
!  Unit Name   : regolith_reworking_zone
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Simulate a lunar regolith reworking zone by vertical mixing of several layers in our push-pop system        
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
subroutine regolith_reworking_zone(user,surf,finterval)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_reworking_zone
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   real(DP),intent(in) :: finterval

   ! Regotrack Internals
   integer(I4B) :: i,j,k
   type(regolayertype),pointer :: current
   type(regolayertype) :: newlayer
   ! Mixing 
   real(DP) :: zmix1mm, zmixbig, mixtime
   real(DP) :: zmix
   real(DP) :: z, z0, zmare, ztot
   integer(I4B), parameter :: LUN=7
   real(DP), parameter :: t0 = 5.64099e-07
   real(DP), parameter :: dts = 0.908976 ! slope in depth vs time plot from Arnold (1975)
   character(len=255), parameter :: filename = 'mixtime.dat'
   logical :: exist
   integer(I4B) :: mixfreq, bigratio
   
   ! Executable code
   ! Reading a time information
   inquire(file=filename, exist=exist)
   if (exist) then 
      open(LUN, file=filename, status='old')
      read(LUN,*) mixtime
   else
      write(*,*) filename,' is missing!'
   end if
   close(LUN)

   mixtime = mixtime + finterval * user%interval

   open(LUN, file=filename, status='replace')
   write(LUN,*) mixtime
   close(LUN)

   ! Apply Arnold's distrubed depth formula
   !zmix1mm = t0 * ( 1.0e+06 )**(dts) / 100.0_DP ! Convert unit from centimeters in original formula to meters
   !zmixbig = t0 * ( mixtime )**(dts) / 100.0_DP
   !zmix     = t0 * (finterval * user%interval)**(dts) / 100.0_DP
   zmix = 0.00017_DP
   
   !bigratio = int(mixtime/1.0e+07)
   !if (mod(bigratio,10) == 0 .and. bigratio /= 0) then 
   !   mixfreq = 2 
   !   zmix(1) = zmix1mm
   !   zmix(2) = zmixbig
   !else
   !mixfreq = 1
   !zmix = zmixbig
   !zmix = 0.01_DP
   !zmix = zmix1mm
   !end if

   !write(*,*) mixtime, bigratio, mod(bigratio,10), mixfreq, zmixbig
   !write(*,*) mixtime, finterval, zmix

   do j=1,user%gridsize
      do i=1,user%gridsize
  
         current => surf(i,j)%regolayer
         z = surf(i,j)%regolayer%thickness

         if (z <= zmix) then 
            z0 = 0._DP
            zmare = 0._DP
            ztot = 0._DP

            do

             if (.not. associated(current%next)) exit

             if (z <= zmix) then
                ztot  = ztot  + current%thickness
                zmare = zmare + current%thickness * current%comp
                current => current%next
                z0 = z
                z = z + current%thickness
             else 
                ztot  = ztot  + (zmix - z0)
                zmare = zmare + (zmix - z0) * current%comp
                exit
             end if

            end do

            call regolith_traverse_pop(-1.0_DP * ztot, surf(i,j))
            newlayer%thickness  = ztot
            newlayer%comp       = zmare/ztot
            newlayer%meltfrac   = 0.0_DP
            call regolith_push(surf(i,j),newlayer)

         end if

      end do
   end do

   return
end subroutine regolith_reworking_zone
