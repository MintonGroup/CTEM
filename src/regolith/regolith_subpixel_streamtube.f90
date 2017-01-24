!****f* regolith/regolith_subpixel_streamtube
! Name
!   regolith_subpixel_streamtube -- Calculate a segment of a stream tube under resolution. 
! SYNOPSIS
!   This uses
!   * module_globals
!   * module_regolith
!
!   call regolith_subpixel_streamtube()
!
! DESCRIPTION
!   
!   A stream tube that is under resolution means that the origin and the emerging place of 
!   a stream tube is in the same pixel. Therefore, the total volume of a stream tube can be 
!   calculated with no need of calculation of a segment. Yet, a stream tube within a pixel 
!   may occupy several layers, and the stream tube will contain different components from 
!   different layers at the pixel. This subroutine is to estimate the slices of a stream tube.
!  
! ARGUMENTS
!   Input
!   * user      -- The user-defined variables from the input file 
!   * surfi     -- A given pixel from surface grid
!   * deltar    -- The size of a stream tube
!   * ri        -- A point of a stream tube's projection on surface grid
!   * rip1      -- A point of a stream tube's projection on surface grid
!   * eradi     -- The inner radial distance of a stream tube
!
!   Output
!   * newlayer  -- 
! 
!***

!**********************************************************************************************************************************
!
!  Unit Name   : regolith_subpixel_streamtube
!  Unit Type   : subroutine
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Subpixel stream tube's approximation 
!  
!
!  Input
!    Arguments :   
!
!  Output
!    Arguments : surf      ::  surface 
!           
! 
!  Notes       :  
!
!**********************************************************************************************************************************
subroutine regolith_subpixel_streamtube(user,surfi,deltar,ri,rip1,eradi,vseg,newlayer,rm)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_subpixel_streamtube
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(inout) :: surfi
   real(DP),intent(in)            :: deltar,ri,rip1,eradi,vseg
   type(regodatatype),intent(inout) :: newlayer
   real(DP),intent(in) :: rm

   ! Traversing a linked list 
   real(DP),parameter :: a = 0.936457 
   real(DP),parameter :: b = 1.12368
   type(regolisttype),pointer :: current
   real(DP) :: z,zmax,zstart,zend,rlefti,rleftf,rrighti,rrightf,rc,vsgly,x

   ! The depth that a stream tube dips 
   zmax = rip1/4.0_DP

   current => surfi%regolayer
   z = surfi%regolayer%regodata%thickness
   zstart = 0.0_DP
   zend = z 

   ! Two cases: subpixel is inside the first layer, and its volume is simply the landing ejecta blanket.
   if (zend>=zmax) then 

      newlayer%thickness = vseg
      newlayer%comp      = vseg * surfi%regolayer%regodata%comp

   else
   ! The subpixel stream tube may dip deeper than the first layer. And we use the line segments approximation, but 
   ! the layer now intersects with both sides of a stream tube, so two intersection points between the layer and the 
   ! stream tube must be calculated. 
   !                                                                                  
   !-----|--------------------------------------------------------|---------------*------------------------------*--- > 0
   !  *  |                                                        |           *                              *
   !  ** |                                                        |    *                              *
   !-----*--------------------------------------------------------*----------------------------*--------------------- > z
   !    *|  *                                              *      |
   !     |     *                                     *            |                    *
   !     |  *      *                         *                    |            *
   !     |                *     *     *                           |    * 
   !     |     *                                                  |
   !     |                                                   *    |
   !     |         *                                   *          |
   !     |               *                     *                  |
   !     |                      *     *                           |
   !^    ^                      ^                                 ^              ^
   !0.0 rleftf                  rc                              rrighti        eradi

     rlefti  = 0.0_DP
     rleftf  = 0.0_DP
     rrighti = eradi
     rrightf = eradi
     rc      = rip1 * sqrt(3.0) / 4.0
     newlayer%thickness = 0._DP
     newlayer%comp      = 0._DP

     do 

      ! It should hit the bottom layer before it exits, I think. 
      if (.not. associated(current%next)) exit

      if (zend<zmax) then
 
         rleftf  = regolith_quadratic_func(zend,rip1,rlefti,rc,rlefti)  
         rrighti = regolith_quadratic_func(zend,rip1,rc,rrightf,rrightf)
         vsgly   = 0.25 * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b/eradi * rleftf) + tan(b/eradi * rrightf)&
                 - tan(b/eradi * rlefti) - tan(b/eradi * rrighti)) - abs(b/eradi * rleftf + b/eradi * rrightf &
                 - b/eradi * rlefti - b/eradi * rrighti)) 
       
         newlayer%comp      = newlayer%comp + vsgly * current%regodata%comp
         newlayer%thickness = newlayer%thickness + vsgly

         current => current%next
         z = z + current%regodata%thickness 
         zstart = zend
         zend = z
         rlefti = rleftf
         rrightf = rrighti
      ! final part of a stream tube 
      else 
         vsgly   = 0.25 * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b) - b)) 
         newlayer%comp   = newlayer%comp + (vsgly - newlayer%thickness) * current%regodata%comp
         newlayer%thickness  = vsgly 
         exit
      end if

     end do

   end if

   call regolith_streamtube_head(user,surfi,deltar,newlayer,eradi,rm)

   return
end subroutine regolith_subpixel_streamtube
