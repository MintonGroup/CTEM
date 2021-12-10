#gnuplot 5.0 script


#Common ranges on all of our SFD plots
xmin = 1e-3
xmax = 2.0

ymin = 1e-2
ymax = 1e3

#Common SFD definitions
nemp(r) = 0.0084 * r**(-2)
ngsat(r) = 0.385 * r**(-2)

#This makes labels parallel to the power laws
labx = 10**((log10(xmax) - log10(xmin))*0.98 + log10(xmin))
angscl = (log10(ymax) - log10(ymin)) / (log10(xmax) - log10(xmin))
ang(slp) = atan2(slp,angscl) * 180 / pi


#Define custom color scheme that is colorblind safe
custcolor01 = 'black'
custcolor02 = '#2166ac'
custcolor03 = '#4393c3'
custcolor04 = '#01665e'
custcolor05 = '#8c510a'
custcolor06 = 'black'
custcolor07 = '#8f8686'
custcolor08 = '#343131'
custcolor09 = 'gray90'
custcolor10 = 'gray50'

#Define custom colorblind safe gradient palette
#set palette defined (0 '#014636',1 '#016c59',2 '#02818a',3 '#3690c0',4 '#67a9cf',5 '#a6bddb',6 '#d0d1e6',7 '#ece2f0', 8 '#fff7fb')

set palette defined (0 '#081d58',1 '#253494',2 '#225ea8',3 '#1d91c0',4 '#41b6c4',5 '#7fcdbb',6 '#c7e9b4',7 '#edf8b1',8 '#ffffd9')

#Line styles
set style line 1 lt 1 lc rgb custcolor01 lw 6 pt 6 ps 1
set style line 2 lt 1 lc rgb custcolor02 lw 5 dt (2,1) pt 4 ps 2
set style line 22 lt 1 lc rgb custcolor02 lw 5 dt (1,1) pt 4 ps 2
set style line 3 lt 1 lc rgb custcolor03 lw 5 dt (1,1) pt 4 ps 2
set style line 4 lt 1 lc rgb custcolor04 lw 5 dt (8,3,1,3) pt 4 ps 2
set style line 5 lt 1 lc rgb custcolor05 lw 5 dt (8,3) pt 4 ps 2
set style line 6 lt 1 lc rgb custcolor06 lw 3 dt 1 pt 6 ps 1
set style line 7 lt 1 lc rgb custcolor07 lw 4 
set style line 8 lt 1 lc rgb custcolor08 lw 4 
set style line 9 lt 1 lc rgb custcolor09 lw 3
set style line 10 lt 1 lc rgb custcolor10 lw 3
