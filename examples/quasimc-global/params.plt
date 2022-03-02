#gnuplot 5.0 script
#parameters for the particular run

fe = 5
Kd1 = 0.0315151258095091
psi1 = 2.0
psi2 = 1.25

runname = sprintf('variable-psi-fe%d-%.1f-%g',fe,psi1,psi2)
outname = sprintf('frames/NPF-global-'.runname.'-equilmovie-%06d.'.ftype,ii)

set output outname

eqfile = 'NPFextrap-equilibrium-'.runname.'.dat'
titletext = sprintf("Variable degradation function exponent for f_e = %g\n{/Symbol y}_{small} = %g     {/Symbol y}_{big} = %4.2f",fe,psi1,psi2)
