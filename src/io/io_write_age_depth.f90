!**********************************************************************************************************************************
!
!  Unit Name   : io_write_age_depth
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Writes new files for terrain grids for age information
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
subroutine io_write_age_depth(user,surf,n_size,icrater,ncrat,age,age_2_depth)
   use module_globals
   use module_io, EXCEPT_THIS_ONE => io_write_age_depth
   implicit none

   ! Arguments
   type(usertype),intent(in)                :: user
   type(surftype),dimension(:,:),intent(in) :: surf
   integer(I2B),intent(in)                  :: n_size
   integer(I8B),intent(in)                  :: icrater
   integer(I8B),intent(in)                  :: ncrat
   real(DP),intent(in)                      :: age
   real(SP),dimension(:,:,:),intent(inout)  :: age_2_depth

   ! Regotrack Internals
   integer(I4B) :: i,j,k
   real(SP),dimension(:),allocatable        :: age_prev
   real(SP),dimension(:,:),allocatable      :: agedepthtot
   real(DP)                                 :: depth, depth_prev, recyclratio
   real(DP),dimension(5),parameter          :: sdepth = [0.3, 3.0, 7.0, 8.0, 9.0] ![0.1, 0.5, 0.8, 1.0, 3.0] ![0.3, 2.0, 4.0, 6.0, 10.0]
   type(regolisttype),pointer               :: current => null()
   integer(I4B)                             :: n 

   ! Executable code
   allocate(age_prev(n_size))
   allocate(agedepthtot(n_size,5))

   agedepthtot = 0.0_SP
   
   do k = 1, 5
      do j = 1,user%gridsize
         do i = 1,user%gridsize

            depth       = 0.0
            depth_prev  = 0.0
            current => surf(i,j)%regolayer
            age_prev(:) = current%regodata%age(:) 
            do 
             if (depth > sdepth(k)) then
                recyclratio      = (depth - sdepth(k)) / (depth - depth_prev)
                agedepthtot(:,k) = agedepthtot(:,k) - recyclratio * age_prev(:)
                depth            = sdepth(k)
                exit
             end if
             depth_prev          = depth
             depth               = depth + current%regodata%thickness
             agedepthtot(:,k)    = agedepthtot(:,k) + current%regodata%age(:)
             age_prev(:)         = current%regodata%age(:)
             current => current%next
            end do

         end do
      end do  
   end do
 
   n = max( icrater / ncrat, 1)
   do k = 1, 5
      age_2_depth(n,1,k)          = age
      age_2_depth(n,2:n_size+1,k) = agedepthtot(:,k)
   end do
 
   deallocate(age_prev,agedepthtot)

   return
end subroutine io_write_age_depth
