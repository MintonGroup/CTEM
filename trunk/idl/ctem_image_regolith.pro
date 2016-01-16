pro ctem_image_regolith,ncount,gridsize,pix,regolith,regolith_image
; Generates the regolith depth map image and saves it as a jpeg output image in the 'rego' directory
; outputs dregolith, which may be scaled to use as a console image
Compile_Opt DEFINT32
thisDevice = !D.Name
Set_Plot, 'Z'
Erase
Device, Set_Resolution=[gridsize,gridsize],Set_Pixel_Depth=24, Decomposed=0 
loadct, 39
TVLCT, red, green, blue, /GET

minref = pix * 1.0d-4
regolith_scaled = dblarr(gridsize,gridsize)
maxreg = max(regolith)
minreg = min(regolith)
if minreg lt minref then minreg = minref
if maxreg lt minref then maxreg = minref + 1.0d30
regolith_scaled = regolith > minreg
regolith_scaled = 254.0d0 * ((alog(regolith_scaled) - alog(minreg)) / (alog(maxreg) - alog(minreg)))

; save regolith display
tv, regolith_scaled, 0, 0, xsize=gridsize, ysize=gridsize, /device
regolith_image = TVRD(True=1)
Set_Plot, thisDevice

if (file_test('rego',/DIRECTORY) eq 0) then begin
	file_mkdir,'rego'
endif
fnum = string(ncount,format='(I6.6)')
fname = 'rego/rego' + fnum + '.jpg'
write_jpeg, fname, regolith_image, true=1, quality=100

end
