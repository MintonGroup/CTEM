!****f* regolith/regolith_streamtube_volume_func
! Name
!   regolith_streamtube_volume_func -- Calculate the volume of a segment of a stream tube 
! SYNOPSIS
!   This uses
!   * module_globals
!   * module_regolith
!
!   volume = regolith_streamtube_volume_func()
!
! DESCRIPTION
!   
!   This function is used very frequently if "doregotrack" is turned
!   on. As a a result, we separated it into an independent function.
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

!**********************************************************************************************************
!
!  Unit Name   : regolith_streamtube_volume_func
!  Unit Type   : function
!  Project     : CTEM
!  Language    : Fortran 2003
!
!  Description : Calculate the volume of a segment of stream tube
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
!***********************************************************************************************************
function regolith_streamtube_volume_func(eradi,ri,rip1,deltar) result(vol)
   use module_globals
   use module_regolith, EXCEPT_THIS_ONE => regolith_streamtube_volume_func
   implicit none

   real(DP), intent(in) :: eradi, ri, rip1, deltar
   real(DP)             :: vol
   real(DP), parameter  :: a = 0.936457
   real(DP), parameter  :: b = 1.12368

   vol = 0.25 * PI * (deltar**2) * (a**2) * eradi / b * abs( &
         abs(tan(b/eradi * rip1) - (b/eradi) * rip1) - &
         abs(tan(b/eradi * ri)   - (b/eradi) * ri))

end function regolith_streamtube_volume_func
