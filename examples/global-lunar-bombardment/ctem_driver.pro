pro ctem_driver

;-------------------------------------------------------------------------
;  Jim Richardson, Arecibo Observatory
;  David Minton, Purdue University Dept. of Earth, Atmospheric, & Planetary Sciences
;  July 2014
;  Cratered Terrain Evolution Model display module
;
;  Inputs are read in from the ctem.in file
;
;-------------------------------------------------------------------------
;------------- Initial Setup ----------------
;-------------------------------------------------------------------------
Compile_Opt DEFINT32
!EXCEPT=2

; ----------- input file -----------------------
infilename = 'ctem.in'
DATFILE='ctem.dat'

; ---------- reading input files ----------
seedarr = lon64arr(100)
seedn = 1
totalimpacts = long64(0)
ncount = long64(0)
curyear = 0.d0
restart = "F" 
fracdone = 1.0d0
masstot = 0.d0

ctem_io_read_input,infilename,interval,numintervals,gridsize,pix,seed,numlayers,sfdfile,impfile,maxcrat,ph1,shadedmaxhdefault,shadedminhdefault,shadedminh,shadedmaxh,restart,runtype,popupconsole,saveshaded,saverego,savepres,savetruelist

seedarr(0) = seed
area = (gridsize * pix)^2

;read input data
pnum = file_lines(impfile)
production = dblarr(2,pnum)
productionfunction = dblarr(2,pnum)
openr,LUN,impfile,/GET_LUN
readf,LUN,productionfunction
free_lun,LUN

;create impactor production population
production(0,*) = productionfunction(0,*)
production(1,*) = productionfunction(1,*)*area*interval

;write out corrected production population
openw,1,sfdfile
printf,1,production
close,1
free_lun,1

;set up cratering surface grid and display-only grid
surface_dem = dblarr(gridsize,gridsize)
regolith = dblarr(gridsize,gridsize)

;set up temporary distribution bins
distl = 1
pdistl = 1
odist = dblarr(6,distl)
tdist = dblarr(6,distl)
pdist = dblarr(6,pdistl)
pdisttotal = dblarr(6,pdistl)

datformat = "(I17,1X,I12,1X,E19.12,1X,A1,1X,F9.6,1X,E19.12)"


if strmatch(restart,'F',/fold_case) then begin ; Start with a clean slate
	print, 'Starting a new run'
	curyear = 0.0d0
	totalimpacts = 0
	masstot = 0.d0
	fracdone = 1.0d0

	if strmatch(runtype,'statistical',/fold_case) then begin
		ncount = 1
		openw,LUN,DATFILE,/GET_LUN
		printf,LUN,totalimpacts,ncount,curyear,restart,fracdone,masstot,format=datformat
		for n=0,seedn-1 do begin
			printf,LUN,seedarr(n),format='(I12)'
		endfor
		free_lun,LUN
	endif else begin
		ncount = 0
	endelse

	surface_dem(*,*) = 0.0d0
	regolith(*,*) = 0.0d0

	file_delete, 'tdistribution.dat',/allow_nonexistent

endif else begin ; continue an old run
	print, 'Continuing a previous run'

	ctem_io_read_old,gridsize,surface_dem,regolith,odist,tdist,pdist,mass

	;read in constants file
	openr,LUN,DATFILE,/GET_LUN
	readf,LUN,totalimpacts,ncount,curyear,restart,fracdone,masstot,format=datformat
	seedn = 0
	while ~ eof(LUN) do begin 
		readf,LUN,iseed
		seedarr(seedn) = long64(iseed)
	   seedn=seedn+1
	endwhile	

endelse

openw,FRACDONEFILE,'fracdone.dat',/GET_LUN
openw,REGODEPTHFILE,'regolithdepth.dat',/GET_LUN

;-------------------------------------------------------------------------
; ---------- begin loops ----------
;-------------------------------------------------------------------------
print, 'Beginning loops'

;define number of loop iterations and begin
;numintervals=ceil(endyear/interval)-ncount
while (ncount le numintervals) do begin

	
	; ---------- creating crater population ----------
	if (ncount gt 0) then begin
	
		fnum = string(ncount,format='(I6.6)')
		if (file_test('misc',/DIRECTORY) eq 0) then begin
			file_mkdir,'misc'
		endif
		; save a copy of the ctem.dat file 
		fname = 'misc/ctem_' + fnum + '.dat'
		file_copy, 'ctem.dat', fname, /OVERWRITE

		print, ncount, '  Calling FORTRAN routine'
		;call fortran program to create & count craters
		spawn, './CTEM',/noshell

		; ---------- reading FORTRAN output ----------
		print, ncount, '  Reading FORTRAN output'
		ctem_io_read_old,gridsize,surface_dem,regolith,odist,tdist,pdist,mass

		;read in constants file
		openr,LUN,DATFILE,/GET_LUN
		readf,LUN,totalimpacts,ncount,curyear,restart,fracdone,masstot,format=datformat
		seedn = 0
		while ~ eof(LUN) do begin 
			readf,LUN,iseed
			seedarr(seedn) = long64(iseed)
			seedn = seedn + 1
		endwhile	
		free_lun,LUN

		curyear = curyear + fracdone * interval
		masstot = masstot + mass
		printf,FRACDONEFILE,fracdone,curyear
		flush,FRACDONEFILE

		printf,REGODEPTHFILE,curyear,mean(regolith),max(regolith),min(regolith)
		flush,REGODEPTHFILE

		;save a copy of the binned observed crater distribution
		if (file_test('dist',/DIRECTORY) eq 0) then begin
			file_mkdir,'dist'
		endif
		fname = 'dist/odist_' + fnum + '.dat'
		file_copy, 'odistribution.dat', fname, /OVERWRITE

		; save a copy of the cumulative observed crater distribution 
		fname = 'dist/ocum_' + fnum + '.dat'
		file_copy, 'ocumulative.dat', fname, /OVERWRITE

		;save a copy of the binned true distribution
		fname = 'dist/tdist_' + fnum + '.dat'
		file_copy, 'tdistribution.dat', fname, /OVERWRITE

		;save a copy of the binned idealized production function
		fname = 'dist/pdist_' + fnum + '.dat'
		file_copy, 'pdistribution.dat', fname, /OVERWRITE

		; save a copy of the cumulative true crater distribution if the user requests it
		if strmatch(savetruelist,'T',/fold_case) then begin
			fname = 'dist/tcum_' + fnum + '.dat'
			file_copy, 'tcumulative.dat', fname, /OVERWRITE
		endif

		; save a copy of the impacted mass 
		fname = 'misc/mass_' + fnum + '.dat'
		file_copy, 'impactmass.dat', fname, /OVERWRITE


	endif

	; Get the accumulated production function
	pdisttotal = pdist 
	pdisttotal(3:5,*) = pdist(3:5,*) * curyear / interval

	; ---------- displaying results ----------
	print, ncount, '  Displaying results'

	ctem_image_dem,ncount,gridsize,pix,surface_dem,surface_dem_image
	if strmatch(saverego,'T',/fold_case) then ctem_image_regolith,ncount,gridsize,pix,regolith,regolith_image
	if strmatch(saveshaded,'T',/fold_case) then begin
		if (shadedminhdefault eq 1) then shadedminh = min(surface_dem)
		if (shadedmaxhdefault eq 1) then shadedmaxh = max(surface_dem)
		ctem_image_shaded_relief,ncount,gridsize,pix,surface_dem,surface_dem,shadedminh,shadedmaxh,'shaded',shaded_image
	endif
	if strmatch(savepres,'T',/fold_case) then ctem_image_presentation,ncount,gridsize,pix,curyear,odist,pdisttotal,tdist,ph1,surface_dem_image
	ctem_window_display,ncount,totalimpacts,gridsize,pix,curyear,masstot,odist,pdisttotal,tdist,ph1,surface_dem,regolith,surface_dem_image,popupconsole

	ncount = ncount + 1

	;write out the current data file
	if (strmatch(runtype,'statistical',/fold_case)) || (ncount eq 1) then begin 
		restart = 'F'
		curyear = 0.0d0
		totalimpacts = 0
		masstot = 0.d0
		file_delete, 'tdistribution.dat',/allow_nonexistent
	endif else begin
		restart = 'T'
	endelse

	openw,LUN,DATFILE,/GET_LUN
	printf,LUN,totalimpacts,ncount,curyear,restart,fracdone,masstot,format=datformat
	for n=0,seedn-1 do begin
		printf,LUN,seedarr(n),format='(I12)'
	endfor
	free_lun,LUN
endwhile
free_lun,FRACDONEFILE
free_lun,REGODEPTHFILE

end
