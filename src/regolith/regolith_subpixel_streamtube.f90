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
subroutine regolith_subpixel_streamtube(user,surfi,deltar,ri,rip1,eradi,newlayer,vmare,totseb,&
                                        age_collector,xmints,xsfints,vol,mixedregodata)
   use module_globals 
   use module_regolith, EXCEPT_THIS_ONE => regolith_subpixel_streamtube
   implicit none

   ! Arguments
   type(usertype),intent(in) :: user
   type(surftype),intent(in) :: surfi
   real(DP),intent(in)       :: deltar,ri,rip1,eradi
   type(regodatatype),intent(inout)    :: newlayer, mixedregodata
   real(DP),intent(out)                :: vmare,totseb
   real(SP),dimension(:),intent(inout) :: age_collector
   real(DP),intent(in)                 :: xmints
   real(DP),intent(in)                 :: xsfints
   real(DP),intent(inout)              :: vol

   ! Internal variables

   ! Traversing a linked list 
   real(DP),parameter :: a = 0.936457 
   real(DP),parameter :: b = 1.12368
   !type(regolisttype),pointer :: current
   type(regodatatype),dimension(:),allocatable :: current
   real(DP) :: z,zmax,zstart,zend,rlefti,rleftf,rrighti,rrightf,rc,vsgly,vsgly1,vsgly2,x,mvl,mvr
   integer(I4B) :: N

   ! Stream tube's distance from the edge of a melt zone 
   real(DP) :: zm, recyratio, xmints1, vseg

   ! Parameters for calculating shocked segment of stream tube   
   real(DP) :: x_up_sh, x_low_sh, vsh

   ! The depth that a stream tube dips 
   zmax = rip1/4.0_DP
   
   !current => surfi%regolayer
   allocate(current,source=surfi%regolayer)
   N = size(current)
   z = surfi%regolayer(N)%thickness
   zstart = 0.0_DP
   zend = z 
   vmare = 0._DP
   totseb = 0._DP
   vseg   = 0.0_DP
   vsgly1  = 0.0_DP
   vsgly2  = 0.0_DP

   ! Two cases: subpixel is inside the first layer, and its volume is simply the landing ejecta blanket.
   if (zend>=zmax) then 
      vmare  = newlayer%thickness * user%pix**2 * surfi%regolayer(N)%comp
      totseb = newlayer%thickness * user%pix**2
      if (eradi>xmints) then
         vseg             = regolith_streamtube_volume_func(eradi,xmints,eradi,deltar)
         vsh              = regolith_shock_damage(eradi,deltar,xmints,xsfints,0.0_DP,eradi)
         mixedregodata%totvolume = vseg-vsh
         recyratio        = max(vseg-vsh,0.0 )/ (user%pix**2) / (surfi%regolayer(N)%thickness)
         mixedregodata%meltvolume = surfi%regolayer(N)%meltfrac * vseg * recyratio
         mixedregodata%meltfrac = mixedregodata%meltvolume / mixedregodata%totvolume
         age_collector(:) = age_collector(:) + surfi%regolayer(N)%age(:) * recyratio
         vol              = vol + sum(age_collector(:))
!         write(*,*) '1',eradi, xmints, xsfints, &
!                    vseg/user%pix**2, (vseg-vsh)/user%pix**2, recyratio
      end if
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
   !     |     *                                     *            |                    *       rrightf
   !     |  *      *                         *                    |            *
   !     |                *     *     *                           |    * 
   !     |     *                                                  |
   !     |                                                   *    |
   !     |         *                                   *          |
   !     |               *                     *                  |
   !     |                      *     *                           |
   !^    ^                      ^                                 ^              ^
   !0.0   rleftf                  rc                              rrighti        eradi
   !rlefti

     rlefti  = 0.0_DP
     rleftf  = 0.0_DP
     rrighti = eradi
     rrightf = eradi
     rc      = rip1 * sqrt(3.0) / 4.0

     do 

      ! It should hit the bottom layer before it exits, I think. 
      !if (.not. associated(current%next)) exit

      if (zend<zmax) then
 
         rleftf  = regolith_quadratic_func(zend,rip1,rlefti,rc,rlefti)  
         rrighti = regolith_quadratic_func(zend,rip1,rc,rrightf,rrightf)
         vsgly1  = regolith_streamtube_volume_func(eradi,rlefti,rleftf,deltar) 
         vsgly2  = regolith_streamtube_volume_func(eradi,rrighti,rrightf,deltar)
         vsgly   = vsgly1 + vsgly2
         vmare = vmare + (vsgly * current(N)%comp)
         totseb = totseb + vsgly
          
         ! If this segmentr, intersecting with layer, is beyond "xmints"
         ! (melt-and-streamline intersection point), a streamtube will retain
         ! whatever it is in the layer. However, if this segment is completely
         ! within shock fragmentation zone, then it would not retain. 
         ! Fragmentation subroutine only should place under the following two
         ! if-then conditions for both segments. If a segment is not satisfied
         ! with the following if-then conditions, it means that a segment is
         ! inside the melt zone that is previously calculated in
         ! "regolith_melt_glass" subroutine. Then next step is to examine how a
         ! segment locates in our shock pressure decay zone. The rule of thumb
         ! is a segment remaining the same volume no matter what it experiences
         ! shock pressues. We should only examine how the ends of a streamtube's
         ! segment locate inside our decay zone. In the first if-then condition,
         ! it is to describe the segment closer to emerging location. Its ends
         ! are "rrighti" and "rrightf". The volume is "vsgly2". The "xsh" is the
         ! x-location of a current layer intersecting with shock pressure decay
         ! zone.
         ! 1) rrighti > xsh:           no damage, retaining ratio (recyratio) is 1.0
         ! 2) rrighti < xsh < rrightf: partial damage, retaining ratio is
         !                             affected to be (vsgly2 - v_shocked) / vsgly2.
         ! 3) rrightf < xsh:           complete damage, retaining ratio is 0.0
         !if (rrighti>max(xmints,sqrt(rm**2 - z**2))) then
         if (rrighti > xmints) then
            vsh                = regolith_shock_damage(eradi,deltar,xmints,xsfints,rrighti,rrightf)
            recyratio          = max((vsgly2 -vsh),0.0)/ (user%pix**2) / current(N)%thickness
            age_collector(:)   = age_collector(:) + current(N)%age(:) * recyratio
            vol                = vol + sum(current(N)%age(:)) * recyratio
            mvr                = vsgly2 * current(N)%meltfrac * recyratio
         end if

         if (rlefti > xmints) then
            vsh                = regolith_shock_damage(eradi,deltar,xmints,xsfints,rlefti,rleftf)
            recyratio          = max((vsgly1-vsh),0.0) / (user%pix**2) / current(N)%thickness
            age_collector(:)   = age_collector(:) + current(N)%age(:) * recyratio
            vol                = vol + sum(current(N)%age(:)) * recyratio
            mvl                = vsgly1 * current(N)%meltfrac * recyratio
         end if

         !current => current%next
         mixedregodata%meltvolume = mixedregodata%meltvolume + mvl + mvr
         mixedregodata%totvolume = mixedregodata%totvolume + vsgly
         mixedregodata%meltfrac = mixedregodata%meltvolume / mixedregodata%totvolume
         if (mixedregodata%meltfrac > 1.0_DP) then
            write(*,*) "ERROR! mixedregodata%meltfrac >1! (SUBPIXEL)"
         end if
         N = N - 1
         z = z + current(N)%thickness 
         zstart = zend
         zend = z
         rlefti = rleftf
         rrightf = rrighti
      ! final part of a stream tube 
      else !I think this means it's in the melt zone like Ya-Huei said above.. if so, nothing needed here.
         vsgly   = 0.25 * PI * deltar**2 * a**2 * eradi / b * (abs(tan(b) - b)) 
         vmare   = vmare + (vsgly - totseb) * current(N)%comp
         totseb  = vsgly 
         exit
      end if
     end do
   end if

   call regolith_streamtube_head(user,surfi,deltar,vmare,totseb,age_collector, mixedregodata)

   return
end subroutine regolith_subpixel_streamtube
