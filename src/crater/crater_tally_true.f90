!**********************************************************************************************************************************
!
!  Unit Name   : crater_tally_true
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Tallies the craters and returns the final distributions
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
subroutine crater_tally_true(domain,truelist,ntrue,truedist)
   use module_globals
   use module_io
   use module_util
   use module_crater, EXCEPT_THIS_ONE => crater_tally_true
   implicit none

   ! Arguments
   type(domaintype),intent(in) :: domain
   integer(I4B),intent(in) :: ntrue
   real(DP),dimension(TRUECOLS,ntrue),intent(inout) :: truelist
   real(DP),dimension(:,:),intent(out)  :: truedist

   ! Internal variables
   integer(I4B) :: i,craternum
   integer(I4B),dimension(ntrue) :: ind
   real(DP),dimension(TRUECOLS,ntrue) :: truetmp


   ! Executable code

   ! Reset all the distribution bins
   do i=1,domain%distl
      truedist(1,i) = 1e3_DP*SQRT2**(domain%plo+(i-1))
      truedist(2,i) = 1e3_DP*SQRT2**(domain%plo+i)
      truedist(3,i) = 0.0_DP 
      truedist(4,i) = 0.0_DP
      truedist(5,i) = 0.0_DP
      truedist(6,i) = 0.0_DP
   end do

   ! Bin the true crater distribution
   do craternum = 1,ntrue
      ! Find out what bin this crate belongs in
      i = ceiling(log(truelist(1,craternum)/1e3_DP)/LOGSQRT2) - domain%plo
      if (i < 1) then
         write(*,*) 'fcrat = ',truelist(1,craternum)
         write(*,*) 'smallest bin: ',1e3_DP*SQRT2**(domain%plo)
         write(*,*) 'domain%smallest_crater = ',domain%smallest_crater
         write(*,*) 'domain%subcrater_limit = ',domain%subcrater_limit
         read(*,*)
      end if
      truedist(3,i) = truedist(3,i) + log(truelist(1,craternum)) ! Geometric mean (intermediate step)
      truedist(4,i) = truedist(4,i) + 1 ! Differential number
   end do

   do i = 1,domain%distl
      if (truedist(4,i) > 0._DP) then
         truedist(3,i) = exp(truedist(3,i) / truedist(4,i)) ! Geometric mean (final step)
      else
         truedist(3,i) = sqrt(truedist(1,i) * truedist(2,i))
      end if
      truedist(5,i) = sum(truedist(4,i:domain%distl)) ! Cumulative number
      truedist(6,i) = (truedist(4,i)*truedist(3,i)**3)/(domain%area*(truedist(2,i)-truedist(1,i))) ! R-value
   end do

   ! Replace the list of all craters with a sorted list
   truetmp(:,1:ntrue) = truelist(:,1:ntrue)
   call util_mrgrnk(truelist(1,1:ntrue),ind(1:ntrue))
   do i = 1,ntrue
      truelist(:,i) = truetmp(:,ind(ntrue - i + 1))
   end do
   return
end subroutine crater_tally_true

