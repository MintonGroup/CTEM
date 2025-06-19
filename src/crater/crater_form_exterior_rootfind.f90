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


subroutine crater_form_exterior_rootfind(user,surf,crater,domain,deltaMtot)
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_form_exterior_rootfind
   implicit none

  ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),dimension(:,:),intent(inout) :: surf
   type(cratertype),intent(inout) :: crater
   type(domaintype),intent(inout) :: domain
   real(DP),intent(in) ::   deltaMtot

   ! Internal variables

   logical :: lastloop 
   real(DP) :: deltaMp,rd
   real(DP) :: x1,x2,startrd
   real(DP),parameter :: Tolerance = 1.e-18_DP
   integer(I4B),parameter :: maxIterations = 100
   real(DP) :: valueAtRoot
   real(DP),parameter :: FPP = 1.e-11_DP
   real(DP),parameter :: nearzero = 1.e-20_DP
   integer(I4B) :: niter,error
   real(DP) :: resultat,AA,BB,CC,DD,EE,FA,FB,FC,Tol1,PP,QQ,RR,SS,xm
   integer(I4B) :: i,j, ev,br
   integer(I4B),parameter :: NTRY = 50
   integer(I4B),parameter :: NBRACKET = 10
   real(DP),parameter :: FIRSTFACTOR = 3.2_DP
   real(DP) :: factor
   real(DP) :: f1,f2

   ! Executable code
   if (abs(deltaMtot) < VSMALL) return
   if (deltaMtot > 0._DP) then 
      write(*,*) 'Too much mass already! Skipping exterior raised rim.',crater%fcratpx
      return
   end if
   lastloop = .false.
   factor = FIRSTFACTOR
   startrd = RIMDROP


   ! First bracket the root
   bracket: do br=1,NBRACKET
      factor = factor / 2
      x1 = startrd
      x2 = startrd / factor
      f1 = crater_form_exterior_func(user,surf,crater,domain,x1,deltaMtot,lastloop)
      f2 = crater_form_exterior_func(user,surf,crater,domain,x2,deltaMtot,lastloop)
      do j = 1, NTRY
         if (f1 * f2 < 0._DP) exit bracket
         if (abs(f1) < abs(f2)) then
            x1 = x1 + factor * (x1 - x2)
            f1 = crater_form_exterior_func(user,surf,crater,domain,x1,deltaMtot,lastloop)
         else
            x2 = x2 + factor * (x2 - x1)
            f2 = crater_form_exterior_func(user,surf,crater,domain,x2,deltaMtot,lastloop)
         end if
         if (abs(f1) > huge(f1).or.(abs(f2) > huge(f2))) cycle bracket
         !write(*,*)
         !write(*,*) j,factor,x1,x2,f1,f2,crater%fcrat
         !read(*,*)
      end do
   end do bracket
   !!if (f1 * f2 >= 0._DP) then
   !   write(*,*) crater%fcrat
   !   write(*,*) deltaMtot
   !   read(*,*)
   !end if

   ! Now do a Brent's method to find the root
   error = 0
   AA = x1
   BB = x2
   FA = f1
   FB = f2
   CC = AA; FC = FA; DD = BB - AA; EE = DD
   if (.not.util_rootbracketed(FA,FB)) then 
      error = -1
      resultat = x1
      write(*,*) 'Error in crater_form_exterior_rootfind! Root not bracketed!'
   else 
      FC = FB 
      do i = 1,maxIterations !,while (done == 0.and.i < maxIterations)
         niter = i
         if (.not.util_rootbracketed(FC,FB)) then
            CC = AA; FC = FA; DD = BB - AA; EE = DD
         end if
         if (abs(FC) < abs(FB)) then
            AA = BB; BB = CC; CC = AA
            FA = FB; FB = FC; FC = FA
         end if
         Tol1 = 2 * FPP * abs(BB) + 0.5_DP * Tolerance
         xm = 0.5_DP * (CC-BB)
         if ((abs(xm) <= Tol1).or.(abs(FB) < nearzero)) then
            ! A root has been found
            rd = BB 
            lastloop = .true.
            deltaMp = crater_form_exterior_func(user,surf,crater,domain,rd,deltaMtot,lastloop) 
            exit
         else 
            if ((abs(EE) >= Tol1).and.(abs(FA) > abs(FB))) then
               SS = FB / FA 
               if (abs(AA - CC) < nearzero) then
                  PP = 2._DP * xm * SS 
                  QQ = 1._DP - SS 
               else 
                  QQ = FA / FC 
                  RR = FB / FC 
                  PP = SS * (2._DP * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1._DP)) 
                  QQ = (QQ - 1._DP) * (RR - 1._DP) * (SS - 1._DP) 
               end if
               if (PP > nearzero) QQ = -QQ 
               PP = abs(PP) 
               if ((2._DP*PP)<min(3._DP*xm *QQ-abs(Tol1*QQ),abs(EE*QQ))) then
                  EE = DD   
                  DD = PP / QQ 
               else 
                  DD = xm   
                  EE = DD 
               end if
            else 
               DD = xm 
               EE = DD 
            end if
            AA = BB 
            FA = FB 
            if (abs(DD) > Tol1) then 
              BB = BB + DD 
            else 
               if (xm > 0) then 
                  BB = BB + abs(Tol1)
               else 
                  BB = BB - abs(Tol1)
               end if
            end if
            FB = crater_form_exterior_func(user,surf,crater,domain,BB,deltaMtot,lastloop)  
         end if
      end do
      if (i >= maxIterations) error = -2
   end if
   !write(*,'(I4,4F19.12)') niter,rd,RIMDROP,deltaMp
   !read(*,*)
   return


end subroutine crater_form_exterior_rootfind

