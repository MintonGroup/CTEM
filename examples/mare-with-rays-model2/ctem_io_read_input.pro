pro ctem_io_read_input,infilename,interval,numintervals,gridsize,pix,seed,numlayers,sfdfile,impfile,maxcrat,ph1,shadedmaxhdefault,shadedminhdefault,shadedminh,shadedmaxh,restart,runtype,popupconsole,saveshaded,saverego,savepres,savetruelist
Compile_Opt DEFINT32
print, 'Reading input file'
openr,infile,infilename, /GET_LUN
line=""
comment="!"
interval = 0.d0
numintervals = 0
pix=-1.0d0 
gridsize=-1
seed = 0
maxcrat = 1.0d0
shadedmaxhdefault = 1
shadedminhdefault = 1
shadedminh = 0.d0
shademaxnh = 0.d0


; Set required strings to unset value
notset="-----NOTSET----"
sfdfile = notset
impfile = notset
sfdcompare = notset
restart = notset
runtype = notset
popupconsole = notset
saveshaded = notset
saverego = notset
savepres = notset
savetruelist = notset
while (not EOF(infile)) do begin
	readf,infile,line
	if (~strcmp(line,comment,1)) then begin
		substrings = strsplit(line,' ',/extract)
		if strmatch(substrings(0),'pix',/fold_case) then reads,substrings(1),pix
		if strmatch(substrings(0),'gridsize',/fold_case) then reads,substrings(1),gridsize
		if strmatch(substrings(0),'seed',/fold_case) then reads,substrings(1),seed
		if strmatch(substrings(0),'sfdfile',/fold_case) then reads,substrings(1),sfdfile
		if strmatch(substrings(0),'impfile',/fold_case) then reads,substrings(1),impfile
		if strmatch(substrings(0),'maxcrat',/fold_case) then reads,substrings(1),maxcrat
		if strmatch(substrings(0),'sfdcompare',/fold_case) then reads,substrings(1),sfdcompare
		if strmatch(substrings(0),'interval',/fold_case) then reads,substrings(1),interval
		if strmatch(substrings(0),'numintervals',/fold_case) then reads,substrings(1),numintervals
		if strmatch(substrings(0),'popupconsole',/fold_case) then reads,substrings(1),popupconsole
		if strmatch(substrings(0),'saveshaded',/fold_case) then reads,substrings(1),saveshaded
		if strmatch(substrings(0),'saverego',/fold_case) then reads,substrings(1),saverego
		if strmatch(substrings(0),'savepres',/fold_case) then reads,substrings(1),savepres
		if strmatch(substrings(0),'savetruelist',/fold_case) then reads,substrings(1),savetruelist
		if strmatch(substrings(0),'runtype',/fold_case) then reads,substrings(1),runtype
		if strmatch(substrings(0),'restart',/fold_case) then reads,substrings(1),restart
		if strmatch(substrings(0),'shadedminh',/fold_case) then begin
			reads,substrings(1),shadedminh
			shadedminhdefault = 0
		endif
		if strmatch(substrings(0),'shadedmaxh',/fold_case) then begin
			reads,substrings(1),shadedmaxh
			shadedmaxhdefault = 0
		endif
	end
end
if interval le 0.0d0 then begin
	print,'Invalid value for or missing variable INTERVAL in ' + infilename
	stop
end
if numintervals le 0 then begin
	print,'Invalid value for or missing variable NUMINTERVALS in ' + infilename
	stop
end
if pix le 0.0d0 then begin
	print,'Invalid value for or missing variable PIX in ' + infilename
	stop
end
if gridsize le 0 then begin
	print,'Invalid value for or missing variable GRIDSIZE in ' + infilename
	stop
end
if seed eq 0 then begin
	print,'Invalid value for or missing variable SEED in ' + infilename
	stop
end
if strmatch(sfdfile,notset) then begin
	print,'Invalid value for or missing variable SFDFILE in ' + infilename
	stop
end
if strmatch(impfile,notset) then begin
	print,'Invalid value for or missing variable IMPFILE in ' + infilename
	stop
end
if strmatch(popupconsole,notset) then begin
	print,'Invalid value for or missing variable POPUPCONSOLE in ' + infilename
	stop
end
if strmatch(saveshaded,notset) then begin
	print,'Invalid value for or missing variable SAVESHADED in ' + infilename
	stop
end
if strmatch(saverego,notset) then begin
	print,'Invalid value for or missing variable SAVEREGO in ' + infilename
	stop
end
if strmatch(savepres,notset) then begin
	print,'Invalid value for or missing variable SAVEPRES in ' + infilename
	stop
end
if strmatch(savetruelist,notset) then begin
	print,'Invalid value for or missing variable SAVETRUELIST in ' + infilename
	stop
end
if strmatch(runtype,notset) then begin
	print,'Invalid value for or missing variable RUNTYPE in ' + infilename
	stop
end
if strmatch(restart,notset) then begin
	print,'Invalid value for or missing variable RESTART in ' + infilename
	stop
end

free_lun,infile

ph1 = dblarr(3,1)
if ~strmatch(sfdcompare,notset) then begin
	cnum = file_lines(sfdcompare)
	ph1 = dblarr(3,cnum)
	openr,COMP,sfdcompare,/GET_LUN
	readf,COMP,ph1
	close,COMP
	free_lun,COMP
end
free_lun,infile

end
