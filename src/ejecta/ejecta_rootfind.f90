!**********************************************************************************************************************************
!
!  Unit Name   : ejecta_rootfind
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Brent's method to search for erad as a function of lrad
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
subroutine ejecta_rootfind(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
   use module_globals
   use module_util
   use module_ejecta, EXCEPT_THIS_ONE => ejecta_rootfind
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   type(domaintype),intent(in) :: domain
   real(DP),intent(inout) :: erad
   real(DP),intent(in) :: lrad
   real(DP),intent(out) :: vejsq,ejang
   logical,intent(inout) :: firstrun

   ! Internals
   real(DP) :: x2,starterad
   real(DP),parameter :: Tolerance = 1.e-7_DP
   integer(I4B),parameter :: maxIterations = 100
   real(DP) :: valueAtRoot
   real(DP),parameter :: FPP = 1.e-11_DP
   real(DP),parameter :: nearzero = 1.e-20_DP
   integer(I4B) :: niter,error
   real(DP) :: resultat,AA,BB,CC,DD,EE,FA,FB,FC,Tol1,PP,QQ,RR,SS,xm
   integer(I4B) :: i,j, ev,br
   integer(I4B),parameter :: NTRY = 50
   integer(I4B),parameter :: NBRACKET = 5
   real(DP),parameter :: FIRSTFACTOR = 1.6_DP
   real(DP) :: factor
   real(DP) :: f1,f2
   integer(I4B) :: numev,numbr

   ! Executable code
   factor=FIRSTFACTOR
   numev=0
   numbr=0
   starterad=erad 
   everything: do ev=1,NTRY
      numev=numev+1
      if (ev==NTRY) then
         write(*,*) 'Major failure in ejecta_rootfind! '
         return
      end if
      bracket: do br=1,NBRACKET
         numbr=numbr+1
         erad = starterad
         x2 = starterad/factor

         ! First bracket the root
         f1 = ejecta_blanket_func(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
         f2 = ejecta_blanket_func(user,crater,domain,x2,lrad,vejsq,ejang,firstrun)
         do j=1,NTRY
            if (f1*f2 < 0._DP) exit
            if (abs(f1) < abs(f2)) then
               erad=erad+factor*(erad-x2)
               f1=ejecta_blanket_func(user,crater,domain,erad,lrad,vejsq,ejang,firstrun)
            else
               x2=x2+factor*(x2-erad)
               f2=ejecta_blanket_func(user,crater,domain,x2,lrad,vejsq,ejang,firstrun)
            end if
         end do
         if ((erad<=crater%rad).and.(erad>0._DP)) exit
         factor=0.5_DP*(factor+1._DP)
      end do bracket

      ! Now do a Brent's method to find the root
      error = 0
      AA = erad
      BB = x2
      FA = f1
      FB = f2
      CC = AA; FC = FA; DD = BB - AA; EE = DD
      if (.not.util_rootbracketed(FA,FB)) then 
         error = -1
         resultat = erad
         write(*,*) 'Error in ejecta_rootfind! Root not bracketed!'
      else 
         FC = FB 
         do i=1,maxIterations !,while (done == 0.and.i < maxIterations)
            if (.not.util_rootbracketed(FC,FB)) then
               CC = AA; FC = FA; DD = BB - AA; EE = DD
            end if
            if (abs(FC) < abs(FB)) then
               AA = BB; BB = CC; CC = AA
               FA = FB; FB = FC; FC = FA
            end if
            Tol1 = 2 * FPP * abs(BB) + 0.5_DP * Tolerance
            xm = 0.5_DP * (CC-BB)
            if ((abs(xm) <= Tol1).or.(abs(FA) < nearzero)) then
               ! A root has been found
               resultat = BB 
               valueAtRoot = ejecta_blanket_func(user,crater,domain,resultat,lrad,vejsq,ejang,firstrun)
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
               FB = ejecta_blanket_func(user,crater,domain,BB,lrad,vejsq,ejang,firstrun)
            end if
         end do
         if (i >= maxIterations) error = -2
      end if
      niter = i
      erad = resultat
      if ((erad <= crater%frad).and.(lrad >= crater%frad).and.(erad >= 0._DP)) exit
      factor = 0.5_DP  * (factor + 1.0_DP) ! Failed. Try again with a new factor
   end do everything

   return
end subroutine ejecta_rootfind

