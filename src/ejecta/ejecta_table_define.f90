!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_table_define
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Computes the ejecta blanket look-up table 
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
subroutine ejecta_table_define(user,crater,domain,ejb,ejtble,melt)
   use module_globals
   use module_util
   use module_regolith 
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_table_define
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(inout)    :: domain
   type(ejbtype),dimension(EJBTABSIZE),intent(out) :: ejb
   integer(I4B),intent(out) :: ejtble
   real(DP),intent(out),optional :: melt

   ! Internal variables
   integer(I4B) :: k
   real(DP) :: erad,eradold,thick,vejsq,ejang,lrad
   logical :: firstrun

   ! Regotrack internal variables
   real(DP) :: rmelt,depthb,dimp,vimp

   ! This will be replaced by its own function
   real(DP) :: c0,c1,c2,c3,flrad,rh,fld,r
   rh = crater%rimheight 
   fld = -crater%floordepth 
   flrad = 0.5_DP * crater%floordiam / crater%frad 
   c1 = (fld - rh) / (flrad + flrad**2 / 3._DP - flrad**3 / 6._DP - 7._DP / 6._DP)
   c0 = rh - (7._DP / 6._DP) * c1
   c2 = c1 / 3._DP
   c3 = -c2 / 2._DP
   !^^^^^^^^^^^^^^^

   ! Executable code

   ! Get estimate of size of ejb table
   crater%continuous = RCONT * crater%frad**(EXPCONT)  ! Continuous ejecta distance From Moore (1974) eq. 1
   crater%ejdis = DISEJB * crater%continuous
                                                        ! We go out a factor of 3 to get the discontinuous ejecta thickness 
   domain%ejbres = (log(crater%ejdis) - log(crater%ejrad)) / EJBTABSIZE
   lrad = crater%ejrad !exp(log(crater%rad) !+ domain%ejbres)
   erad = crater%ejrad 
   ejtble = EJBTABSIZE
   firstrun = .true.
   thick = 0._DP
   
   if (present(melt)) then
    
      if (user%testflag) then 
         dimp = user%testimp
         vimp = user%testvel
      else
         dimp = crater%imp
         vimp = crater%impvel
      end if 
   

      call regolith_melt_zone(user,crater,dimp,vimp,rmelt,depthb)
   end if
 
   !write(*,*) '     lrad/Df                   vej                       ebh                       melt & 
   !           fraction              melt thickness'
   do k = 0,EJBTABSIZE
      call ejecta_rootfind(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
      if (k >= 1) then
         !call ejecta_thickness(user,crater,eradold,erad,lrad - domain%ejbres,lrad,thick)
         ejb(k)%lrad = log(lrad)

         ! This will be replaced 
         r = lrad / crater%frad
         if (lrad >= crater%frad) then
            thick = crater%rimheight * r**(-3.0_DP)
         else
            !thick = 0.14_DP * crater%frad**(0.74_DP) / (crater%frad - crater%ejrad) * (lrad - crater%ejrad)
            thick = c0 + c1 * r + c2 * r**2 + c3 * r**3
         end if
         ejb(k)%thick = log(thick) 
         ejb(k)%vesq = vejsq
         ejb(k)%angle = ejang
         ejb(k)%erad = log(erad)
         if (present(melt)) then
            call regolith_melt_fraction(dimp,depthb,erad,eradold,rmelt,melt)
            ejb(k)%meltfrac = melt
         end if
         !write(*,*) lrad/crater%rad,erad/crater%rad,sqrt(vejsq),thick !,melt,thick*melt
         if ((thick <= VSMALL) .or. (abs(eradold - erad) < VSMALL)) then
            ejtble = k
            crater%ejdis = lrad
            exit
         end if
      end if
      lrad = exp(log(lrad) + domain%ejbres)
      eradold = erad
   end do
   !write(*,*) 'A MELT ZONE of ',crater%frad,' meter-sized crater: ',rmelt,'at a rim',ejb(1)%meltfrac
   ! Get pixel space distance
   crater%ejdispx = nint(crater%ejdis / user%pix)

   return
end subroutine ejecta_table_define
