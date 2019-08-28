pro ctem_image_dem,ncount,gridsize,pix,surface_dem,surface_dem_image
; Generates the shaded surface DEM image and saves it as a jpeg output image in the 'surf' directory
; outputs surface_dem_arr, which may be scaled to use as a console image

; This code is for reading in the x-y positions from the cumulative distribution and separating out into layers
; so that the tallied craters can be drawn as circles
;!PATH = Expand_Path('+/home/campus/daminton/coyote/') + ':' + !PATH

Compile_Opt DEFINT32
thisDevice = !D.Name
Set_Plot, 'Z'
Erase
Device, Set_Resolution=[gridsize,gridsize],Set_Pixel_Depth=24, Decomposed=0 
TVLCT, red, green, blue, /GET

sun = 20.d0
radsun = sun * !dtor

surface_dem_arr = dblarr(gridsize,gridsize)

surface_dem_arr=surface_dem-shift(surface_dem,0,1)
surface_dem_arr  = (0.5d0*!dpi) + atan(surface_dem_arr,pix)  ; Get average slope
surface_dem_arr = abs(surface_dem_arr - radsun < 0.5d0*!dpi) ; incident angle
surface_dem_arr = 254.0d0  * cos(surface_dem_arr)            ; shaded relief surface

tv, surface_dem_arr, 0, 0, xsize=gridsize, ysize=gridsize, /device
print,'plotting'

surface_dem_image = TVRD(True=1)
Set_Plot, thisDevice

; save surface display
if (file_test('surf',/DIRECTORY) eq 0) then begin
	file_mkdir,'surf'
endif
fnum = string(ncount,format='(I6.6)')
fname = 'surf/surf' + fnum + '.jpg'

write_jpeg, fname, surface_dem_image, true=1, quality=100

end
