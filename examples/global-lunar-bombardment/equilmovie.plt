#gnuplot 5.0 script
do for [ii=166:166] {
reset


ftype = "eps"
set terminal postscript enhanced eps color size 6.5,3.65625 18
load "params.plt"
load "common.plt"
show output


set multiplot layout 1,2 title titletext font ",18"
set size square

set format x ""
set format y ""

set xrange[000:1999]
set yrange[000:1999]

unset xtics 
unset ytics

sg = 1.0
lg = 3.7

set tmargin 1.0
set bmargin 2.1
set lmargin sg
set rmargin lg

pix = 3.08e3
gridsize = 2000

set obj 1 rect from 1600,30 to 1600+100e4/pix,60 fc rgb "white" front fs solid 1.0 noborder
set label 1 "1000 km" at 1600+50e4/pix,100 tc rgb "white" front center

set style line 1 lt 1 lc rgb "white" lw 0.1

surf = sprintf("surf/surf%06d.jpg",ii)
circ = sprintf("dist/ocum_%06d.dat",ii)

plot surf binary filetype=auto w rgbimage notitle

unset obj 1
unset label 1

xmin = 10.0
xmax = 1000.0

ymin = 1e-8
ymax = 1e-3

set logscale xy

set format x "10^{%L}"
set format y "10^{%L}"

set xtics 10
set ytics 10

set mxtics 10
set mytics 10

set rmargin sg+0.5
set lmargin lg

set xlabel "Crater diameter {/Times-Italic D_c} (km)"
set ylabel "Cumulative number per unit area {/Times-Italic n_{>D_c}} (km^{-2})" offset -1

#set style line 1 lt 1 lc rgb "black" lw 1 pt 6 ps 1
#set style line 2 lt 1 lc rgb "black" lw 3 dt 2 pt 4 ps 2
#set style line 3 lt 1 lc rgb "black" lw 3 dt 3 pt 4 ps 2
#set style line 4 lt 1 lc rgb "black" lw 3 dt 4 pt 4 ps 2
#set style line 5 lt 1 lc rgb "black" lw 3 dt 6 pt 4 ps 2
#set style line 9 lt 1 lc rgb "blue" lw 3

A = (pix * gridsize)**2

pdist = sprintf("dist/pdist_%06d.dat",ii)
ocum = sprintf("dist/ocum_%06d.dat",ii)




set xrange[xmin:xmax]
set yrange[ymin:ymax]

labx = 10**((log10(xmax) - log10(xmin))*0.97 + log10(xmin))
angscl = (log10(ymax) - log10(ymin)) / (log10(xmax) - log10(xmin))

ang(slp) = atan2(slp,angscl) * 180 / pi

set style textbox opaque fillcolor "white" noborder
set label 9  "Geometric saturation" at first labx, first 8.0*ngsat(labx) rotate by ang(-2.0) right boxed front tc rgb custcolor05
set label 10 "Predicted equilibrium" at first labx, first 2*ngsat(labx) rotate by ang(-2.0) right boxed front tc rgb custcolor02

set fit quiet
set fit logfile '/dev/null'
f(x) = a*x+b
fit [log10(20):log10(80)] [*:*] f(x) eqfile u (log10($1*2*1e-3)):(log10($2*1e6))  via a,b 

# get the relation of x- and y-range
dx = log10(xmax)-log10(xmin)
dy = log10(ymax)-log10(ymin)
s1 = dx/dy
# get ratio of axes
# helper function for getting the rotation angle of the labels in degree
deg(x) = x/pi*180.0
r(x) = deg(atan(s1*x))


#set label 11 "Empirical equilibrium" at first labx, first 1.9*nemp(labx) rotate by ang(-2.0) center boxed front

Xval = sprintf("{/Times-Italic X }= %5.3f",ii / 600.0)
set label 20 Xval at screen 0.98, screen 0.92 right

set key spacing 1.3 at screen 0.83, screen 0.30 box width -3 opaque
#print labx

Ahighlands = 20051255.514 / 1.5
AMoon = 3.793e7
AFassett = 2.798701
set style line 17 lt 1 lc rgb custcolor06 lw 1 pt 7 ps 1.00
set style line 27 lt 1 lc rgb custcolor07 lw 1 pt 8 ps 1.00
set style line 37 lt 1 lc rgb custcolor05 lw 1 pt 9 ps 1.00

load "NPF.plt"

if (ii > 0) { 
   plot \
    "NonMareNonSpa-Highlands.csv" u ($1 < 200 ? $1 : NaN):(($0+1)/Ahighlands):(sqrt($0+0.9999999)/Ahighlands) w errorbars ls 27 title "LOLA Highlands",\
    "Neumann2015GRAIL_farside_basins.dat" u 1:(2*($0+1)/AMoon):(2*sqrt($0+0.999999)/AMoon) w errorbars ls 37 title "GRAIL Far side basins",\
    pdist u ($1*1e-3):($5/A*ii*1e6) w lines ls 4 title "CTEM Production",\
    eqfile u ($1*2*1e-3):($2*1e6) w lines ls 2 notitle "Equilibrium",\
    ngsat(x*0.5e3)*1e6 w lines ls 5 notitle,\
    ocum u ($1*1e-3):(($0+1)/A*1e6) w points ls 6 title "CTEM Counts"


} else {
   plot \
    "../NonMareNonSpa-Highlands.csv" u ($1 < 200 ? $1 : NaN):(($0+1)/Ahighlands):(sqrt($0+0.9999999)/Ahighlands) w errorbars ls 27 title "LOLA Highlands",\
    "../Neumann2015GRAIL_farside_basins.dat" u 1:(2*($0+1)/AMoon):(2*sqrt($0+0.999999)/AMoon) w errorbars ls 37 title "GRAIL Far side basins",\
    NaN w lines ls 4 title "CTEM Production",\
    eqfile u ($1*2*1e-3):($2*1e6) w lines ls 2 notitle "Equilibrium",\
    ngsat(x*0.5e3)*1e6 w lines ls 5 notitle,\
    NaN w points ls 6 title "CTEM Counts"
}

     #"" u ($0 == 291 ? $1*2*1e-3 : NaN):(4.0*$2*1e6):(sprintf("Predicted equilibrium")) w labels rotate by ang(-1.8) right boxed tc rgb custcolor02 notitle,\
     #"" u ($0 == 291 ? $1*2*1e-3 : NaN):(4.0*$2*1e6):(sprintf("Predicted equilibrium")) w labels rotate by ang(-1.8) right boxed tc rgb custcolor02 notitle,\

makeanim = sprintf("mogrify -verbose -format png -density 600 -flatten -resize 1920x1080\! ".outname)
system(sprintf(makeanim." ; rm ".outname."; cd .."))

unset multiplot
}
