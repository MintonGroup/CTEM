pro ctem_image_presentation,ncount,gridsize,pix,curyear,odist,pdist,tdist,ph1,map
; Generates the simplified version of the console display for use in presentations
; Plots the R-plot on the left and the image map on the right

; presentation graphics
presresx = 1024
presresy = 0.6*presresx
area = (gridsize * pix)^2

minx = (pix / 3.0d0) *  1d-3 
maxx = 3 * pix * gridsize * 1d-3

;strings for display
dum = string(format='(E10.4)', curyear)
parts = strsplit(dum, 'E', /extract)
fcuryear = strcompress(string(format='(F6.4,"x10!U",i,"!N")', parts[0], parts[1]))
time = 'Time = '
timeunit = ' yr'

if (file_test('presentation',/DIRECTORY) eq 0) then begin
	file_mkdir,'presentation'
endif

thisDevice = !D.Name
Set_Plot, 'Z'
Erase
Device, Set_Resolution=[presresx,presresy],Set_Pixel_Depth=24, Decomposed=0 
loadct, 0
TVLCT, red, green, blue, /GET

;geometric saturation
geom = dblarr(2,2)
geomem = dblarr(2,2)
geomep = dblarr(2,2)
geomel = dblarr(2,2)
geom(0,0) = minx
geom(0,1) = maxx
geom(1,*) = 3.12636d0
geomem(0,0) = minx
geomem(0,1) = maxx
geomem(1,*) = 0.156318d0
geomep(0,0) = minx
geomep(0,1) = maxx
geomep(1,*) = 0.312636d0
geomel(0,0) = minx
geomel(0,1) = maxx
geomel(1,*) = 0.0312636d0

; Remove zeros
nz = where(odist(5,*) ne 0.0,count)
if (count gt 0) then begin
	odistnz = dblarr(6,count)
	odistnz(*,*)= odist(*,nz)
endif else begin
	odistnz = odist
endelse

nz = where(tdist(5,*) ne 0.0,count)
if (count gt 0) then begin
	tdistnz = dblarr(6,count)
	tdistnz(*,*)= tdist(*,nz)
endif else begin
	tdistnz = tdist
endelse

nz = where(pdist(5,*) ne 0.0,count)
if (count gt 0) then begin
	pdistnz = dblarr(6,count)
	pdistnz(*,*)= pdist(*,nz)
endif else begin
	pdistnz = pdist
endelse

; create r-plot array containing exactly 1 crater per bin
area = (gridsize * pix * 1d-3)^2
plo = 1
while (sqrt(2.0d0)^plo gt minx) do begin
	plo = plo - 1
endwhile
phi = plo + 1
while (sqrt(2.0d0)^phi lt maxx) do begin
	phi = phi + 1
endwhile
n = phi - plo
sdist = dblarr(6,n + 1)
p = plo
for i=0,n  do begin
	sdist(0,i) = sqrt(2.d0)^p 
	sdist(1,i) = sqrt(2.d0)^(p+1) 
	sdist(2,i) = sqrt(sdist(0,i) * sdist(1,i))
	sdist(3,i) = 1.0d0
	sdist(5,i) = (sdist(2,i))^3 / (area * (sdist(1,i) - sdist(0,i)))
	p = p + 1
endfor
sdist(4,*) = reverse(total(sdist(3,*),/cumulative),2)


plot, odistnz(2,*)*1.0d-3, odistnz(5,*), line=0, color=255, thick=2.0, $
   TITLE='Crater Distribution R-Plot', charsize=1.2, $
   XTITLE='Crater Diameter (km)', $
   XRANGE=[minx,maxx], /XLOG, XSTYLE=1, $
   YTITLE='R Value', $
   YRANGE=[5.0e-4,5.0e0], /YLOG, YSTYLE=1, $
   /device, pos = [0.12,0.12,0.495,0.495]*presresx

Dfac = sqrt(sqrt(2.0d0))

;display observed crater distribution
oplot, tdistnz(2,*)*1.0d-3, tdistnz(5,*), line=0, color=255, thick=1.0
;oplot, geom(0,*), geom(1,*), line=2, color=255, thick=1.0
oplot, geomem(0,*), geomem(1,*), line=1, color=255, thick=1.0
oplot, geomep(0,*), geomep(1,*), line=1, color=255, thick=1.0
;oplot, geomel(0,*), geomel(1,*), line=1, color=255, thick=1.0
oplot, pdistnz(2,*)*1.0d-3, pdistnz(5,*), line=3, color=255, thick=1.0
oplot, odistnz(2,*)*1.0d-3, odistnz(5,*), line=1, color=255, thick=1.0
oplot, sdist(0,*), sdist(5,*), line=1, color=255, thick=1.0
oplot, ph1(0,*)*Dfac*1.0d-3, ph1(1,*), psym=1, color=255, thick=1.0
;oplot, ph2(0,*), ph2(1,*), psym=4, color=255, thick=1.0
;oplot, ph3(0,*), ph3(1,*), psym=5, color=255, thick=1.0

;draw box around the main surface display & fill
displaysize =  floor(0.45*presresx) ; size of displayed surface
surfxpos = floor(0.520*presresx)
surfypos = floor(0.12*presresy)
xbox = [surfxpos - 1,surfxpos - 1,surfxpos + displaysize,surfxpos + displaysize,surfxpos - 1]
ybox = [surfypos - 1,surfypos + displaysize,surfypos + displaysize,surfypos - 1,surfypos - 1]
plotS, xbox, ybox, /device, color=255

mapscaled = congrid(map,3,displaysize,displaysize,/interp)

;mapscaled=congrid(map,displaysize,displaysize) ; scale image
tv, mapscaled, surfxpos, surfypos, xsize=displaysize, ysize=displaysize, true=1, /device

xyouts, 0.310*presresx, 0.04*presresy, time + fcuryear + timeunit, color=255, charsize=2*presresy/720., /device

snapshot = TVRD(True=1)
Set_Plot, thisDevice
fnum = string(ncount,format='(I6.6)')
fname = 'presentation/presentation' + fnum + '.jpg'
Write_JPEG, fname, snapshot, True=1, Quality=90


end
