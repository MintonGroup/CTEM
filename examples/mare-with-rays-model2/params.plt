#gnuplot 5.0 script
#parameters for the particular run

fe = 10
Kd1 = 0.00312669649143281
psi1 = 2.0
psi2 = 1.4

runname = sprintf('raymodel2')
outname = sprintf('frames/mare-with-rays-'.runname.'-equilmovie-%06d.'.ftype,ii)

set output outname

titletext = sprintf("Minton et al. (2019) Ray Model 2")
