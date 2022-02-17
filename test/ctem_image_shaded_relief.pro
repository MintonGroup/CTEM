pro ctem_image_shaded_relief,ncount,gridsize,pix,surface_dem,height_vals,minh,maxh,dirname,shaded_image
; Generates the shaded depth map image and saves it as a jpeg output image in the 'rego' directory
; outputs dshaded, which may be scaled to use as a console image
; Uses the array height_vals for the color
Compile_Opt DEFINT32
thisDevice = !D.Name
Set_Plot, 'Z'
Erase
Device, Set_Resolution=[gridsize,gridsize],Set_Pixel_Depth=24, Decomposed=0 
loadct, 61;72; 17 ;33
TVLCT, red, green, blue, /GET

light=[[1,1,1],[0,0,0],[-1,-1,-1]]

; convolution of the dem with the 3x3 matrix

shaded=bytscl(convol(float(surface_dem),float(light)))
if max(shaded) eq 0 then shaded=255

; scale the dem to the height
if (minh gt maxh) then begin
	tmp = minh
	minh = maxh
	maxh = tmp
endif

demscaled = ((height_vals - minh)  > 0.0) < (maxh-minh)
if ((maxh - minh) eq 0.0) then begin
	demscaled = demscaled * 0.0
endif else begin
	demscaled = demscaled/(maxh-minh)*255.0
endelse

tv, demscaled, 0, 0, xsize=gridsize, ysize=gridsize, /device
shaded_image = TVRD(True=1)
shadedscl =float(shaded)/255.0
shaded_imagearr = dblarr(3,gridsize,gridsize)
shaded_imagearr(0,*,*) = float(shaded_image(0,*,*)) * shadedscl(*,*)
shaded_imagearr(1,*,*) = float(shaded_image(1,*,*)) * shadedscl(*,*)
shaded_imagearr(2,*,*) = float(shaded_image(2,*,*)) * shadedscl(*,*)
shaded_image=round(shaded_imagearr)
Set_Plot, thisDevice

if (file_test(dirname,/DIRECTORY) eq 0) then begin
	file_mkdir,dirname
endif

fnum = string(ncount,format='(I6.6)')
fname = dirname + '/' + dirname + fnum + '.jpg'
write_jpeg,fname,shaded_image,true=1,quality=90

end
