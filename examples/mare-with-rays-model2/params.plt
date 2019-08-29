#gnuplot 5.0 script
#parameters for the particular run

pix = 1.0
gridsize = 1000

fe = 10.0
Kd1 = 0.00312669649143281
psi = 2.0

eta = 3.2
gamma = 2.0
Kv1 = 0.17
neq1 = 0.0084

fetext = "{/Times-Roman ".(fe < 10 ? gprintf("%.1f",fe) : gprintf("%.0f",fe))."}"
psitext = "{/Times-Roman ".gprintf("%.1f",psi)."}"
etatext = "{/Times-Roman ".gprintf("%.1f",eta)."}"
Kd1text = "{/Times-Roman ".(Kd1 > 0.1 ? gprintf("%.3f",Kd1) : gprintf("%4.3t",Kd1)."{/Symbol \264}10^{".gprintf("%T",Kd1)."}" )."}"

set output sprintf("frames/ray-fe%02.0f-narrow-equilmovie-%06d.".ftype,fe,ii)

titletext = "Crater size-dependent degradation (Ray model 2): {/Times-Italic K_{d,1} }{/Times-Roman = }".Kd1text.";  {/Times-Italic f_e }{/Times-Roman = }".fetext.";  {/Symbol-Oblique y }{/Times-Roman = }".psitext."
