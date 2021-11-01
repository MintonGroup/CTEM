#gnuplot 5.0 script
reset

set terminal postscript enhanced eps color size 5.0,5.0 16
set output 'counttest.eps'


set format x ""
set format y ""
set xrange[0:1999]
set yrange[0:1999]

set size ratio 1.0

unset xtics
unset ytics

pix = 3.08e3

#set label 1 "A" at graph 0.05,graph 0.95 front center tc rgb "white" font ",24"

set style line 1 lt 1 lc rgb "white" lw 1
set style line 2 lt 1 lc rgb "blue" lw 0.1
set style line 3 lt 1 lc rgb "red" lw 0.1
set style line 4 lt 1 lc rgb "green" lw 0.1
set style line 5 lt 1 lc rgb "cyan" lw 0.1


plot "surf/surf000166.jpg" binary filetype=auto w rgbimage notitle,\
     "dist/ocum_000166.dat" u ($4 > 0.00 ? $2/pix - 1 : 1/0):($3/pix - 1):(0.5*$1/pix) w circles ls 1 notitle
