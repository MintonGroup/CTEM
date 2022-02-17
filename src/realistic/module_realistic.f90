module module_realistic
implicit none
public
save

interface
   function realistic_turbulence(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
   use module_globals
   implicit none
   real(DP), intent(in) ::  x, y, noise_height, freq, pers
   integer(I4B) :: num_octaves
   real(DP),dimension(3,num_octaves),intent(in) :: anchor
   real(DP) :: noise
   end function realistic_turbulence
end interface

interface
   function realistic_arrayinput(x, y, num_octaves, Sarr, Aarr, anchor) result(noise)
   use module_globals
   implicit none
   real(DP), intent(in) ::  x, y
   integer(I4B) :: num_octaves
   real(DP), dimension(num_octaves),intent(in) :: Sarr, Aarr
   real(DP), dimension(3,num_octaves),intent(in) :: anchor
   real(DP) :: noise
   end function realistic_arrayinput
end interface

interface
   function realistic_billowedNoise(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
   use module_globals
   implicit none
   real(DP), intent(in) ::  x, y, noise_height, freq, pers
   integer(I4B) :: num_octaves
   real(DP),dimension(3,num_octaves),intent(in) :: anchor
   real(DP) :: noise
   end function realistic_billowedNoise
end interface

interface
   function realistic_plawNoise(x, y, noise_height, freq, pers, slope, num_octaves, anchor) result(noise)
   use module_globals
   implicit none
   real(DP), intent(in) ::  x, y, noise_height, freq, pers, slope
   integer(I4B) :: num_octaves
   real(DP),dimension(3,num_octaves),intent(in) :: anchor
   real(DP) :: noise
   end function realistic_plawNoise
end interface

interface
   function realistic_ridgedNoise(x, y, noise_height, freq, pers, num_octaves, anchor) result(noise)
   use module_globals
   implicit none
   real(DP), intent(in) ::  x, y, noise_height, freq, pers
   integer(I4B) :: num_octaves
   real(DP),dimension(3,num_octaves),intent(in) :: anchor
   real(DP) :: noise
   end function realistic_ridgedNoise
end interface

interface
   function realistic_swissTurbulence(x, y, lacunarity, gain, warp, num_octaves, anchor) result(noise)
   use module_globals
   implicit none
   real(DP), intent(in) ::  x, y, lacunarity, gain, warp
   integer(I4B) :: num_octaves
   real(DP),dimension(3,num_octaves) :: anchor
   real(DP) :: noise
   end function realistic_swissTurbulence
end interface

interface
   function realistic_jordanTurbulence(x,y,lacunarity,gain1,gain,warp0,warp,damp0,damp,damp_scale,&
           num_octaves, anchor) result(noise)
   use module_globals
   implicit none
   real(DP),intent(in) :: x, y, lacunarity, gain1, gain, warp0, warp, damp0, damp, damp_scale
   integer(I4B),intent(in) :: num_octaves
   real(DP),dimension(3,num_octaves),intent(in) :: anchor
   real(DP) :: noise
   end function realistic_jordanTurbulence
end interface

interface
   subroutine realistic_perlin_noise(xx,yy,noise,dx,dy)
   use module_globals
   implicit none
   real(DP),intent(in) :: xx,yy
   real(DP),intent(out) :: noise
   real(DP),intent(out),optional ::  dx, dy
   end subroutine realistic_perlin_noise
end interface


interface
      subroutine realistic_ejecta_texture(user,surf,crater,deltaMtot,inc,ejecta_dem)
      use module_globals
      implicit none
      type(usertype),intent(in) :: user
      type(surftype),dimension(:,:),intent(inout) :: surf
      type(cratertype),intent(inout) :: crater
      real(DP),intent(in) :: deltaMtot
      integer(I4B),intent(in) :: inc
      real(DP),dimension(-inc:inc,-inc:inc),intent(inout) :: ejecta_dem
      end subroutine realistic_ejecta_texture
end interface

end module module_realistic

