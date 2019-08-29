pro ctem_io_read_old,gridsize,surface_dem,regolith,odist,tdist,pdist,mass
Compile_Opt DEFINT32
openr,LUN,'surface_ejc.dat',/GET_LUN
readu,LUN,regolith
free_lun,LUN

openr,LUN,'surface_dem.dat',/GET_LUN
readu,LUN,surface_dem
free_lun,LUN

distl = file_lines('odistribution.dat') - 1
pdistl = file_lines('pdistribution.dat') - 1

;read in observed cumulative distribution
odist = dblarr(6,distl)
line = "temp"
openr,LUN,'odistribution.dat',/GET_LUN
readf,LUN,line
readf,LUN,odist
close,LUN
free_lun,LUN

;read in true cumulative distribution
tdist = dblarr(6,distl)
openr,LUN,'tdistribution.dat',/GET_LUN
readf,LUN,line
readf,LUN,tdist
close,LUN
free_lun,LUN

;read in production function
pdist = dblarr(6,pdistl)
openr,LUN,'pdistribution.dat',/GET_LUN
readf,LUN,line
readf,LUN,pdist
close,LUN
free_lun,LUN

;read in accumulated mass from file
openr,LUN,'impactmass.dat',/GET_LUN
readf,LUN,mass
close,LUN
free_lun,LUN

end
