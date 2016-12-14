pro ctem_window_display,ncount,totalimpacts,gridsize,pix,curyear,masstot,odist,pdist,tdist,ph1,surface_dem,regolith,map,popupconsole
; Produces the console window display. The array 'map' is used to generate the image on the right-hand side
Compile_Opt DEFINT32

; console output graphics resolution
conresx = 1280
;conresy = 0.758*conresx
profsize = 0.00
conresy = (0.60 + profsize) * conresx 

minx = (pix / 3.0d0) *  1d-3 
maxx = 3 * pix * gridsize * 1d-3
charfac = 1.6 ; character size multiplier

if strmatch(popupconsole,'T',/fold_case) then begin
	device, decomposed=0, retain=2
	thisDevice = !D.NAME
	SET_PLOT, 'Z'
	TVLCT, red, green, blue, /GET
   SET_PLOT, thisDevice
	!P.MULTI = [0, 1, 2]
	Window, 0, xsize=conresx, ysize=conresy
endif else begin
	SET_PLOT, 'Z', /COPY
	device, set_resolution=[conresx,conresy],decomposed=0, Z_Buffer=0
	TVLCT, red, green, blue, /GET
	!P.MULTI = [0, 1, 2]
	erase
endelse

;set up first color scale
loadct, 0

;strings for display
time = 'Time = '
timeunit = ' yr'
timp = 'Total impacts (dotdash) = '
tcrt = 'Total craters (line) = '
ocrt = 'Countable craters (bold) = '
epwr = 'Mean regolith depth = '
mtot = 'Impacted mass = '
c1lab = 'Min. Elevation'
c2lab = 'Mean Elevation'
c3lab = 'Max. Elevation'

maxelev = max(surface_dem)
minelev = min(surface_dem)
medelev = mean(surface_dem)
c1 = string(minelev, format = '(F9.1)') + ' m'
c2 = string(medelev, format = '(F9.1)') + ' m'
c3 = string(maxelev, format = '(F9.1)') + ' m'
tslp = mean(regolith)

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
sq2 = sqrt(2.0d)
while (sq2^plo gt minx) do begin
	plo = plo - 1
endwhile
phi = plo + 1
while (sq2^phi lt maxx) do begin
	phi = phi + 1
endwhile
n = phi - plo
sdist = dblarr(2, n + 1)
p = plo
for i=0, n  do begin
  sdist(0,i) = sq2^p
  sdist(1,i) = sq2^(2.0 * p + 1.5)/ (area * (sq2 - 1))
  p = p + 1
endfor

plot, odistnz(2,*)*1.0d-3, odistnz(5,*), line=0, color=255, thick=2.0, $
   TITLE='Crater Distribution R-Plot', charsize=charfac*conresy/720, $
   XTITLE='Crater Diameter (km)', $
   XRANGE=[minx,maxx], /XLOG, XSTYLE=1, $
   YTITLE='R Value', $
   YRANGE=[5.0e-4,5.0e0], /YLOG, YSTYLE=1, $
   /device, pos = [0.085*conresx,(0.25 + profsize)*conresy,0.44*conresx,0.85*conresy]

Dfac = sqrt(sqrt(2.0d0))

;display observed crater distribution
oplot, tdistnz(2,*)*1.0d-3, tdistnz(5,*), line=0, color=255, thick=1.0
oplot, geom(0,*), geom(1,*), line=2, color=255, thick=1.0
oplot, geomem(0,*), geomem(1,*), line=1, color=255, thick=1.0
oplot, geomep(0,*), geomep(1,*), line=1, color=255, thick=1.0
oplot, geomel(0,*), geomel(1,*), line=1, color=255, thick=1.0
oplot, pdistnz(2,*)*1.0d-3, pdistnz(5,*), line=3, color=255, thick=1.0
oplot, odistnz(2,*)*1.0d-3, odistnz(5,*), line=1, color=255, thick=1.0
oplot, sdist(0,*), sdist(1,*), line=1, color=255, thick=1.0
oplot, ph1(0,*)*Dfac*1.0d-3, ph1(1,*), psym=1, color=255, thick=1.0
;oplot, ph2(0,*), ph2(1,*), psym=4, color=255, thick=1.0
;oplot, ph3(0,*), ph3(1,*), psym=5, color=255, thick=1.0

; set up profile display
;surfpro = dblarr(gridsize)
;surfhpro = dblarr(gridsize)
;surfrpro = dblarr(gridsize)
;surfpos = dindgen(gridsize)
;surfpos = (-0.5d0*(gridsize-1) + surfpos) * pix
;surfpro(*) = surface_dem(*,gridsize/2-1)
;surfhpro(*) = regolith(*,gridsize/2-1)
;surfrpro(*) = surfpro(*) - surfhpro(*)

; set up profile
;maxelevp = max(surfpro)
;minelevp = min(surfpro)
;medelevp = mean(surfpro)
;vertrange = 1.5 * abs(maxelevp-minelevp)
;vertadd = 0.5d0 * (vertrange - abs(maxelevp-minelevp))
;horzrange = vertrange * 5.86666666667d0
;if (horzrange gt (pix*gridsize)) then horzrange = pix*gridsize
;plot, surfpos(*)/1.0d3, surfpro(*)/1.0d3, line=0, color=255, thick=1.0, $
;  TITLE='Matrix Cross-Section', charsize=1.0, $
;  XTITLE='Horizontal Position (km)', $
;  XRANGE=[(-1.0d0*horzrange/2.0d3),(horzrange/2.0d3)], XSTYLE=1, $
;  YTITLE='Elevation (km)', $
;  YRANGE=[((minelevp - vertadd)/1.0d3),((maxelevp + vertadd)/1.0d3)], YSTYLE=1, $
;  /device, pos = [0.063*conresx,0.049*conresy,0.99*conresx,0.26*conresy]
;oplot, surfpos(*)*1.0d-3, surfrpro(*)*1.0d-3, line=0, color=255, thick=1.0

;display text
dum = string(format='(E10.4)', curyear)
parts = strsplit(dum, 'E', /extract)
fcuryear = strcompress(string(format='(F6.4,"x10!U",i,"!N")', parts[0], parts[1])) + ' yr'
ftimp = string(totalimpacts, format = '(I13)')
ftcnt = string(tdist(4,0), format = '(I13)')
focnt = string(odist(4,0), format = '(I13)')
ftslp = string(tslp, format = '(F13.3)') + ' m'


dum = string(format='(E10.4)',masstot)
parts = strsplit(dum, 'E', /extract)
fmtot = strcompress(string(format='(F6.4,"x10!U",i,"!N")', parts[0], parts[1])) + ' kg'

; Display information text below the R-plot
texttop = 0.15 + profsize
textspc = 0.035
textleft = 0.05*conresx
textbreak = 0.25*conresx
xyouts, textleft, (texttop - 0*textspc)*conresy, time, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft, (texttop - 1*textspc)*conresy, timp, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft, (texttop - 2*textspc)*conresy, tcrt, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft, (texttop - 3*textspc)*conresy, ocrt, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft, (texttop - 4*textspc)*conresy, mtot, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft + textbreak, (texttop - 0*textspc)*conresy, fcuryear, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft + textbreak, (texttop - 1*textspc)*conresy, ftimp, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft + textbreak, (texttop - 2*textspc)*conresy, ftcnt, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft + textbreak, (texttop - 3*textspc)*conresy, focnt, color=255, charsize=charfac*conresy/720., /device
xyouts, textleft + textbreak, (texttop - 4*textspc)*conresy, fmtot, color=255, charsize=charfac*conresy/720., /device

; Display information text above the R-plot
xyouts, 0.0368*conresx, 0.972*conresy, c1lab, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.195*conresx, 0.972*conresy, c2lab, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.353*conresx, 0.972*conresy, c3lab, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.0368*conresx, 0.944*conresy, c1, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.195*conresx, 0.944*conresy, c2, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.353*conresx, 0.944*conresy, c3, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.100*conresx, 0.917*conresy, epwr, color=255, charsize=charfac*conresy/720., /device
xyouts, 0.278*conresx, 0.917*conresy, ftslp, color=255, charsize=charfac*conresy/720., /device

;print to screen
print, ncount, '  time = ', curyear
;print, ncount, '  tot impacts = ', tdist
;print, ncount, '  tot craters = ', totalimpacts
;print, ncount, '  obs craters = ', ocrcnt
;print, ncount, '  regolith = ', tslp
;print, ncount, '  prod R = ', mean(prval(1,*))

;draw box around the main surface display & fill
displaysize =  floor(0.53*conresx) ; size of displayed surface
surfxpos = floor(0.462*conresx)
surfypos = floor((0.05 + profsize)*conresy)
xbox = [surfxpos - 1,surfxpos - 1,surfxpos + displaysize,surfxpos + displaysize,surfxpos - 1]
ybox = [surfypos - 1,surfypos + displaysize,surfypos + displaysize,surfypos - 1,surfypos - 1]
plotS, xbox, ybox, /device, color=255

; for display
mapscaled = congrid(map,3,displaysize,displaysize,/interp)

tv, mapscaled, surfxpos, surfypos, xsize=displaysize, ysize=displaysize, true=1, /device

;  ---------- saving window ----------

; save console window
fnum = string(ncount,format='(I6.6)')
; print screen to JPEG file
if (file_test('console',/DIRECTORY) eq 0) then begin
	file_mkdir,'console'
endif
print, ncount, '  saving window'
fname = 'console/console' + fnum + '.jpg'
dwin = tvrd(true=1)
write_jpeg, fname, dwin, true=1, quality=90

end
