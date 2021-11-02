!**********************************************************************************************************************************
!
!  Unit Name   : crater_profile
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Function that defines the basic profile of a crater, not including the ejecta blanket
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
function crater_profile(user,crater,r) result(h)
   use module_globals
   use module_crater, EXCEPT_THIS_ONE => crater_profile
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater
   real(DP),intent(in) :: r

   ! Result variable
   real(DP) :: h

   ! Internal variables
   real(DP) :: flrad,c0,c1,c2,c3
   real(DP),parameter :: A = 4._DP / 11._DP
   real(DP),parameter :: B = -32._DP / 187._DP

   ! Executable code
   flrad = crater%floordiam / crater%fcrat
   
   ! Use polynomial crater profile similar to that of Fassett et al. (2014), but the parameters are set by the crater dimensions
   c1 = (-crater%floordepth - crater%rimheight) / (flrad - 1._DP + A * (flrad**2 - 1._DP) + B * (flrad**3 - 1._DP))
   c0 = crater%rimheight - c1 * (1._DP + A + B)
   c2 = A * c1
   c3 = B * c1

   if (r < flrad) then
      h = -crater%floordepth 
   else if (r > 1.0_DP) then
      h = crater%rimheight * (r**(-RIMDROP))
   else
      h = c0 + c1 * r + c2 * r**2 + c3 * r**3 
   end if

   return
end function crater_profile


function crater_profile_find_r_inner_wall(user,crater) result(r_inner_wall)
   !Uses Brent's method to solve for the point in the crater profile where the inner wall of the crater profile meets the original surface
   use module_globals
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_profile_find_r_inner_wall
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(cratertype),intent(in) :: crater

   ! Result variable
   real(DP) :: r_inner_wall

   integer(I4B),parameter :: maxIterations=100
   real(DP) :: valueAtRoot
 
   real(DP),parameter :: Tolerance = 1e-8_DP
   real(DP),parameter :: FPP = 1.e-11_DP
   real(DP),parameter :: nearzero = 1.e-20_DP
 
   integer :: niter

   real(DP) :: resultat,AA,BB,CC,DD,EE,FA,FB,FC,Tol1,PP,QQ,RR,SS,xm,x1,x2
   integer :: i, done

   i = 0
   done = 0
  
   x1 = 0.0_DP
   x2 = 1.0_DP

   AA = x1
   BB = x2
   FA = crater_profile(user,crater,AA) 
   FB = crater_profile(user,crater,BB) 

   FC = FB;
   do while (done.eq.0.and.i < maxIterations)
      if (.not.util_rootbracketed(FC,FB)) then
         CC = AA; FC = FA; DD = BB - AA; EE = DD
      endif
      if (abs(FC) < abs(FB)) then
         AA = BB; BB = CC; CC = AA
         FA = FB; FB = FC; FC = FA
      endif
      Tol1 = 2 * FPP * abs(BB) + 0.5_DP * Tolerance
      xm = 0.5_DP * (CC-BB)
      if ((abs(xm) <= Tol1).or.(abs(FA) < nearzero)) then
        ! A root has been found
        resultat = BB;
        done = 1
        valueAtRoot = crater_profile(user,crater,resultat) 
      else 
        if ((abs(EE) >= Tol1).and.(abs(FA) > abs(FB))) then
          SS = FB/ FA;
          if (abs(AA - CC) < nearzero) then
            PP = 2 * xm * SS;
            QQ = 1._DP - SS;
          else 
            QQ = FA/FC;
            RR = FB /FC;
            PP = SS * (2 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1._DP));
            QQ = (QQ - 1._DP) * (RR - 1._DP) * (SS - 1._DP);
          endif
          if (PP > nearzero) QQ = -QQ;
          PP = abs(PP);
          if ((2 * PP) < min(3 * xm * QQ - abs(Tol1 * QQ),abs( EE * QQ))) then
            EE = DD;  DD = PP/QQ;
          else 
            DD = xm;   EE = DD;
          endif
        else 
          DD = xm;
          EE = DD;
        endif
        AA = BB;
        FA = FB;
        if (abs(DD) > Tol1) then 
          BB = BB + DD;
        else 
          if (xm > 0._DP) then 
            BB = BB + abs(Tol1)
          else 
            BB = BB - abs(Tol1)
          endif
        endif
        FB = crater_profile(user,crater,BB)
        i=i+1
      endif
   end do
   if (i >= maxIterations) then
      write(*,*) 'Too many iterations in profile solver.'
      r_inner_wall = 1._DP
      return
   end if
   niter = i
   r_inner_wall = resultat
   return
end function crater_profile_find_r_inner_wall

