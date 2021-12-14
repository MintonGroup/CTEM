#gnuplot 5.0 script
# The Neukum production function 
# Ivanov, Neukum, and Hartmann (2001) SSR v. 96 pp. 55-86
# N1 is number of 1 km craters and larger per km^2 vs. time (T)
# N is number of craters per Gy vs. diameter (D) (Neukum et al. (2001), page 60: per km^2?)
# T is Gy
# D is km
# scale factor in time T is N1(T)/N(1.0)

N1lunar(T) = 5.44e-14 * (exp(6.93*T) - 1) + 8.38e-4*T

# Lunar crater SFD
aL00 = -3.0876
aL01 = -3.557528
aL02 =  0.781027
aL03 =  1.021521
aL04 = -0.156012
aL05 = -0.444058
aL06 =  0.019977
aL07 =  0.086850
aL08 = -0.005874
aL09 = -0.006809
aL10 =  8.25e-4
aL11 =  5.54e-5

Nlunar(D) = D < 0.01 ? 1/0 : D > 1000 ? 1/0 : 10**( \
                  aL00 \
					 +	aL01 * log10(D)**1 \
                + aL02 * log10(D)**2 \
                + aL03 * log10(D)**3 \
                + aL04 * log10(D)**4 \
                + aL05 * log10(D)**5 \
                + aL06 * log10(D)**6 \
                + aL07 * log10(D)**7 \
                + aL08 * log10(D)**8 \
                + aL09 * log10(D)**9 \
                + aL10 * log10(D)**10 \
                + aL11 * log10(D)**11 )

#Projectile SFD

aP00 =  0.0
aP01 = -1.375458
aP02 =  1.272521e-1
aP03 = -1.282166
aP04 = -3.074558e-1 
aP05 =  4.149280e-1
aP06 =  1.910668e-1
aP07 = -4.260980e-2
aP08 = -3.976305e-2
aP09 = -3.180179e-3
aP10 =  2.799369e-3
aP11 =  6.892223e-4
aP12 =  2.614385e-6
aP13 = -1.416178e-5
aP14 = -1.191124e-6

Rproj(D) = D < 1e-4 ? 1/0 : D > 300 ? 1/0 : 10**( \
                  aP00 \
					 +	aP01 * log10(D)**1 \
                + aP02 * log10(D)**2 \
                + aP03 * log10(D)**3 \
                + aP04 * log10(D)**4 \
                + aP05 * log10(D)**5 \
                + aP06 * log10(D)**6 \
                + aP07 * log10(D)**7 \
                + aP08 * log10(D)**8 \
                + aP09 * log10(D)**9 \
                + aP10 * log10(D)**10 \
                + aP11 * log10(D)**11 \
                + aP12 * log10(D)**12 \
                + aP13 * log10(D)**13 \
                + aP14 * log10(D)**14 )


# Mars crater SFD
# Ivanov (2001) SSR v. 96 pp. 87-104
aM00 = -3.384
aM01 = -3.197
aM02 = 1.257
aM03 = 0.7915
aM04 = -0.4861
aM05 = -0.3630
aM06 = 0.1016
aM07 = 6.756e-2
aM08 = -1.181e-2
aM09 = -4.753e-3
aM10 = 6.233e-4
aM11 = 5.805e-5

N1mars(T) = 2.68e-14 * (exp(6.93*T) - 1) + 4.13e-4*T

Nmars(D) = D < 0.01 ? 1/0 : D > 1000 ? 1/0 : 10**( \
                  aM00 \
					 +	aM01 * log10(D)**1 \
                + aM02 * log10(D)**2 \
                + aM03 * log10(D)**3 \
                + aM04 * log10(D)**4 \
                + aM05 * log10(D)**5 \
                + aM06 * log10(D)**6 \
                + aM07 * log10(D)**7 \
                + aM08 * log10(D)**8 \
                + aM09 * log10(D)**9 \
                + aM10 * log10(D)**10 \
                + aM11 * log10(D)**11 )
