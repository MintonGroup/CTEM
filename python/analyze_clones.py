import numpy as np
import matplotlib.pyplot as plt
import ctem
import csv
import timeit

def find_apollo_coords(gridsize,pix):
    
    a14xm = -529810
    a14ym = -110559
    a15xm = 110165
    a15ym = 792410
    a16xm = 470012
    a16ym = -272122
    a17xm = 933110
    a17ym = 612289
    
    a14xo = a14xm / pix
    a14yo = a14ym / pix
    a15xo = a15xm / pix
    a15yo = a15ym / pix
    a16xo = a16xm / pix
    a16yo = a16ym / pix
    a17xo = a17xm / pix
    a17yo = a17ym / pix
    
    center = int(gridsize/2)
    
    a14xp = center + int(a14xo)
    a14yp = center + int(a14yo)
    a15xp = center + int(a15xo)
    a15yp = center + int(a15yo)
    a16xp = center + int(a16xo)
    a16yp = center + int(a16yo)
    a17xp = center + int(a17xo)
    a17yp = center + int(a17yo)

    ax = [a14xp,a15xp,a16xp,a17xp]
    ay = [a14yp,a15yp,a16yp,a17yp]

    return ax,ay

def traverse(x,y,depth,rego):
    """finds number of layers to reach depth"""
    n = 1
    nm = 0
    while nm < depth:
        i = -n
        layerdepth = rego[x,y][i]
        nm = nm + layerdepth
        n += 1
    n = n-1
    partial = nm - depth
    thickness = layerdepth - partial
    return n, thickness

def find_melt_at_pixel(x,y,depth,rego,melt,meltdist):
    nlayers, bothickness = traverse(y,x,depth,rego)
    meltarray = np.zeros((nlayers, nlist))
    nm = np.zeros((nlayers, nlist))
    rt = np.zeros(nlayers)
    rt[0] = bothickness
    if nlayers > 1:
        rt[1:] = rego[y,x][-nlayers+1:]
    for i in range(1,nlayers+1):
        if i == 1:
            meltarray[i-1] = meltdist[y,x][-nlist:]
        else:
            meltarray[i-1] = meltdist[y,x][-(i * nlist):-(nlist * (i-1))]
        nm[i-1] = meltarray[i-1] * rt[-i]
    meltarray[nlayers-1] = meltarray[nlayers-1] * (rt[0]/rego[y,x][-nlayers])
    totmelt = np.zeros(nlist)
    for k in range(nlist):
        totmelt[k] = np.sum(meltarray[:,k])
    return totmelt, np.sum(rt)

def aggregate(x,y,depth,n,rego,melt,meltdist,filter=False,apollo=False):
    if n == 0:
        meltarray, rtlist = find_melt_at_pixel(x,y,depth,rego,melt,meltdist)
    elif n < 0:
        print("number of pixels must be >=0")
        return
    else:
        meltarray = []  # might be easiest for this to be a list of arrays?
        rtlist = []
        k = 0
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                meltarray.append(np.zeros(nlist))
                m, rt = find_melt_at_pixel(x+i,y+j,depth,rego,melt,meltdist)
                meltarray[k] = m
                rtlist.append(rt)
                k += 1

    meltarray = np.asarray(meltarray)
    totmelt = np.sum(meltarray)

    norm = np.zeros(nlist)
    totm = np.zeros(nlist)
    for l in range(nlist):
        if n != 0:
            totm[l] = np.sum(meltarray[:,l])
            norm[l] = totm[l] / totmelt
        else:
            totm[l] = meltarray[l]
            norm[l] = totm[l] / totmelt
    mf = totmelt / (np.sum(rtlist) * pix * pix)


    return norm, totmelt, mf
#def analyze_apollo_sites(path,nclones,gridsize,pix):
def analyze_apollo_sites(path,nclones):
	"""Analyzes the Apollo sites for nclones in directories with a suffix in %04g format"""

	apollo_x,apollo_y = find_apollo_coords(gridsize,pix)
	a14normsV = []
	a14meltsV = []
	a14totmeltsV = []
	a14totmfsV = []
	a15normsV = []
	a15meltsV = []
	a15totmeltsV = []
	a15totmfsV = []
	a16normsV = []
	a16meltsV = []
	a16totmeltsV = []
	a16totmfsV = []
	a17normsV = []
	a17meltsV = []
	a17totmeltsV = []
	a17totmfsV = []

	for i in range(0,nclones):
		newpath = path + '%04i' %(i)
		rego = ctem.util.read_linked_list_binary(newpath + '/surface_rego.dat', gridsize)
		melt = ctem.util.read_linked_list_binary(newpath + '/surface_melt.dat', gridsize)
		stack = ctem.util.read_unformatted_binary(newpath + '/surface_stacknum.dat', gridsize, kind='I4B')
		meltdist = ctem.util.read_linked_list_binary(newpath + '/surface_meltdist.dat', gridsize, kind='DP')
		a14normV, a14totmeltV, a14totmfV = aggregate(apollo_x[0],apollo_y[0],depth,aggregation_pixels,rego,melt,meltdist,filter=False)
		a15normV, a15totmeltV, a15totmfV = aggregate(apollo_x[1],apollo_y[1],depth,aggregation_pixels,rego,melt,meltdist,filter=False)
		a16normV, a16totmeltV, a16totmfV = aggregate(apollo_x[2],apollo_y[2],depth,aggregation_pixels,rego,melt,meltdist,filter=False)
		a17normV, a17totmeltV, a17totmfV = aggregate(apollo_x[3],apollo_y[3],depth,aggregation_pixels,rego,melt,meltdist,filter=False)
		a14normsV.append(a14normV)
		a14meltsV.append(a14normV*a14totmeltV)
		a14totmeltsV.append(a14totmeltV)
		a14totmfsV.append(a14totmfV)
		a15normsV.append(a15normV)
		a15meltsV.append(a15normV*a15totmeltV)
		a15totmeltsV.append(a15totmeltV)
		a15totmfsV.append(a15totmfV)
		a16normsV.append(a16normV)
		a16meltsV.append(a16normV*a16totmeltV)
		a16totmeltsV.append(a16totmeltV)
		a16totmfsV.append(a16totmfV)
		a17normsV.append(a17normV)
		a17meltsV.append(a17normV*a17totmeltV)
		a17totmeltsV.append(a17totmeltV)
		a17totmfsV.append(a17totmfV)
		
		

	#save each array
	with open("a14normsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a14normsV)
	with open("a14meltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a14meltsV)
	with open("a14totmeltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a14totmeltsV)
	with open("a14mfsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a14totmfsV)
	with open("a15normsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a15normsV)
	with open("a15meltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a15meltsV)
	with open("a15totmeltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a15totmeltsV)
	with open("a15mfsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a15totmfsV)
	with open("a16normsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a16normsV)
	with open("a16meltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a16meltsV)
	with open("a16totmeltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a16totmeltsV)
	with open("a16mfsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a16totmfsV)
	with open("a17normsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a17normsV)
	with open("a17meltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerows(a17meltsV)
	with open("a17totmeltsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a17totmeltsV)
	with open("a17mfsV.csv", "w") as f:
		wr = csv.writer(f)
		wr.writerow(a17totmfsV)

	return


#######

if __name__ == '__main__':
	start = timeit.default_timer()
	path = '/scratch/bell/blevins2/clone500/ctem_500_clone'
	gridsize = 500
	pix = 12.32e3
	depth = 5
	nlist = 117
	craters = ["South Pole-Aitken", "Fecunditatis", "Australe North", "TOPO-13", "Medii", "Bartels-Voskresenskiy", "Serenitatis North", "Lamont", "Vaporum", "Fowler-Charlier", "Mutus-Vlacq", "Asperitatis", "Aestuum", "Copernicus-H", "Balmer-Kapetyn", "Cruger-Sirsalis", "Dirichlet-Jackson", "Coloumb-Sarton", "Schiller-Zucchius", "Smythii", "Amundsen-Ganswindt", "Nubium", "Rupes Recta", "Deslandres", "Lorentz", "Pozcobutt", "Fitzgerald-Jackson", "Poincare", "Wegener-Winlock", "Orientale Southwest", "TOPO-22", "Galois", "Schickard", "Fermi", "Ingenii", "Gagarin", "Keeler West", "Birkhoff", "Harkhebi", "Campbell", "Serenitatis", "Freundlich-Sharanov", "Landau", "Von Karman M", "Leibnitz", "Apollo", "Milne", "Grimaldi", "Pasteur", "Nectaris", "Mendel-Rydberg", "Moscoviense North", "Moscoviense", "Korolev", "Mendeleev", "Planck", "Hertzsprung", "Humorum", "Crisium", "Humboldtianum", "Sikorsky-Rittenhaus", "Bel'Kovich", "Bailly", "d'Alembert", "Oppenheimer", "Clavius", "Crisium East", "SChwarzschild", "La Condamine", "Maupertuis", "Imbrium", "Schrodinger", "Ukert", "Campanus", "Protagoras", "Calippus", "Babbage A", "Prinz", "Lambert R", "Egede", "Fontenelle", "T. Mayer", "Stadius", "Cassini", "Letronne", "Cruger", "Marius", "Arzachel", "T. Mayer W", "Wallace", "Kepler D", "Orientale", "Lansberg", "Humboldt", "Kies", "Timaeus", "Bonpland D", "Krieger", "Lansberg C", "Briggs", "Archimedes", "Thebit", "Hansteen", "Davy", "Markov", "Iridum", "Seleucus", "Cardanus", "Krafft", "Dalton", "St. George", "Bianchinni", "Herotodus", "Plato", "Damoiseau", "Mairan", "Local"]
	aggregation_pixels = 0
	analyze_apollo_sites(path,50)
	stop = timeit.default_timer()
	execution_time = stop - start
	print("Program Executed in "+str(execution_time)) # It returns time in seconds

