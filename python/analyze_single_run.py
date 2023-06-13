import numpy as np
import matplotlib.pyplot as plt
import ctem
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.image as mpimg


def find_apollo_coords(gridsize, pix):
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

    center = int(gridsize / 2)

    a14xp = center + int(a14xo)
    a14yp = center + int(a14yo)
    a15xp = center + int(a15xo)
    a15yp = center + int(a15yo)
    a16xp = center + int(a16xo)
    a16yp = center + int(a16yo)
    a17xp = center + int(a17xo)
    a17yp = center + int(a17yo)

    axs = [a14xp, a15xp, a16xp, a17xp]
    ays = [a14yp, a15yp, a16yp, a17yp]

    return axs, ays

def traverse(x, y, depth, rego):
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
    if n==0:
        meltarray,rtlist = find_melt_at_pixel(x,y,depth,rego,melt,meltdist)
    elif n<0:
        print("number of pixels must be >=0")
        return
    else:
        meltarray = []
        rtlist = []
        k = 0
        for i in range(-n,n+1):
            for j in range(-n,n+1):
                meltarray.append(np.zeros(nlist))
                m,rt = find_melt_at_pixel(x+i,y+j,depth,rego,melt,meltdist)
                meltarray[k] = m
                rtlist.append(rt)
                k += 1
                
    meltarray = np.asarray(meltarray)
    totmelt = np.sum(meltarray)

    norm = np.zeros(nlist)
    totm = np.zeros(nlist)
    for l in range(nlist):
        if n !=0:
            totm[l] = np.sum(meltarray[:,l])
            norm[l] = totm[l] / totmelt
        else:
            totm[l] = meltarray[l]
            norm[l] = totm[l] / totmelt
    mf = totmelt/(np.sum(rtlist)*pix*pix)
                
        #setup data to plot
    nonzerocraters = []
    nonzeromelt = []
    nonzerocraternames = []
    sigcraters = []
    sigmelt = []
    nonsigmelt = []
    sigcraternames = []
    for l in range(nlist):
        if filter:
            if norm[l] > 1e-4:
                nonzerocraters.append(str(l))
                nonzeromelt.append(norm[l])
                nonzerocraternames.append(craters[l])
                if norm[l] < 0.05:
                    nonsigmelt.append(norm[l])
                else:
                    sigcraters.append(str(l))
                    sigmelt.append(norm[l])
                    sigcraternames.append(craters[l])

        else:
            if norm[l] > 0:
                nonzerocraters.append(str(l))
                nonzeromelt.append(norm[l])
                nonzerocraternames.append(craters[l])
                if norm[l] < 0.05:
                    nonsigmelt.append(norm[l])
                else:
                    sigcraters.append(str(l))
                    sigmelt.append(norm[l])
                    sigcraternames.append(craters[l])

    if len(nonsigmelt) > 0:
        sigcraters.append('%i Others'%len(nonsigmelt))
        sigmelt.append(np.sum(nonsigmelt))
        sigcraternames.append('%i Others'%len(nonsigmelt))

    #plot a pie chart

    f1, (a1, a2) = plt.subplots(1,2)
    f1.suptitle('Pixel: [' + '{:3d}'.format(x) + ',' + '{:3d}'.format(y) + '] (Aggregated to ' + str(n) + ' pixels) (Total melt fraction = ' + '{0:.4f}'.format(mf) + ') (depth = ' + '{0:.2f}'.format(depth) + ' m)')
    a1.pie(sigmelt,labels=sigcraternames,startangle=90,autopct='%1.1f%%',shadow=False,colors = iter(plt.cm.Set3(np.linspace(0,1,len(sigmelt)))))
    a1.axis('equal')
    a2.imshow(img)
    #a2.scatter(x, gridsize-y, color='blue', label='Selected Pixel')
    a2.xaxis.set_ticklabels([])
    a2.yaxis.set_ticklabels([])
    if n == 0:
        a2.scatter(x, gridsize-y, color='blue', label='Selected Region')
    else:
        rect = patches.Rectangle((x,gridsize-y),n,n,linewidth=1,edgecolor='b',facecolor='b',label='Selected Region')
        a2.add_patch(rect)
    if apollo:
        a2.scatter(apollo_x, surfcoords, color='red', label='Apollo 14-17 Landing Sites')
    a2.tick_params(left = False)
    a2.tick_params(bottom = False)
    a2.legend()
    plt.show()

    plt.rcParams['figure.figsize'] = [13, 5]

    return nonzerocraters, nonzerocraternames, nonzeromelt

def map_crater_to_depth(crater,depth,rego,meltdist,melt,apollo=True,log=False):
    if crater > nlist:
        print("Please specify a crater number between 1 and " + str(nlist))
    else:
        cratermelt = np.zeros((gridsize,gridsize))
        cratermf = np.zeros((gridsize,gridsize))
        for x in range(gridsize):
            for y in range(gridsize):
                nlayers, bothickness = traverse(y,x,depth)
                meltarray = np.zeros((nlayers, nlist))
                nm = np.zeros((nlayers, nlist))
                rt = np.zeros(nlayers)
                rt[0] = bothickness
                m = np.zeros(nlayers)
                if nlayers > 1:
                    rt[1:] = rego[y,x][-nlayers+1:]
                    mt = melt[y,x][-nlayers:]
                    mt[0] = mt[0] * (depth/rego[y,x][-nlayers])
                for i in range(1,nlayers+1):
                    if i == 1:
                        meltarray[i-1] = meltdist[y,x][-nlist:]
                    else:
                        meltarray[i-1] = meltdist[y,x][-(i * nlist):-(nlist * (i-1))]
                    nm[i-1] = meltarray[i-1] * rt[-i]
                    meltarray[nlayers-1] = meltarray[nlayers-1] * (rt[0]/rego[y,x][-nlayers])
                tot = np.zeros(nlayers)
                for k in range(nlayers):
                    m[k] = meltarray[k, crater]
                cratermelt[y,x] = np.sum(m)
                cratermf[y,x] = cratermelt[y,x] / (depth * pix * pix)

        f1 = plt.figure(1)
        a1 = f1.add_subplot(111)
        if crater == nlist-1:
            a1.set_title('Local melt to a depth of ' + str(depth) + ' m' )
        else:
            a1.set_title('Melt from crater ' + craters[crater] + ' to a depth of ' + str(depth) + ' m' )
        if not log:
            md = a1.imshow(cratermelt,origin='lower')#,vmin=0.0,vmax=1.0)
        else:
            md = a1.imshow(cratermelt,origin='lower', norm=colors.LogNorm())#vmin=1e-3))#,vmax=1.0))
        if crater == nlist-1:
            cbar = f1.colorbar(md, label='Volume of total melt from local craters')
        else:
            cbar = f1.colorbar(md, label='Volume of total melt from this crater')
        if apollo:
                a1.scatter(apollo_x, apollo_y, color='red', label='Apollo 14-17 Landing Sites')
                plt.legend()
        f2 = plt.figure(2)
        a2 = f2.add_subplot(111)
        if crater == nlist-1:
            a2.set_title('Local melt to a depth of ' + str(depth) + ' m' )
        else:
            a2.set_title('Melt from crater ' + craters[crater] + ' to a depth of ' + str(depth) + ' m' )
        if not log:
            md = a2.imshow(cratermf,origin='lower')#,vmin=0.0,vmax=1.0)
        else:
            md = a2.imshow(cratermf,origin='lower', norm=colors.LogNorm())#vmin=1e-3))#,vmax=1.0))
        if crater == nlist-1:
            cbar = f2.colorbar(md, label='Fraction of total melt from local craters')
        else:
            cbar = f2.colorbar(md, label='Fraction of total melt from this crater')
        if apollo:
                a2.scatter(apollo_x, apollo_y, color='red', label='Apollo 14-17 Landing Sites')
                plt.legend()

        plt.show()

    return

def map_meltfrac_to_depth(depth,apollo=True,log=False):
    mf = np.zeros((gridsize,gridsize))
    for x in range(gridsize):
        for y in range(gridsize):
            nlayers, bothickness = traverse(y,x,depth)
            mfarray = np.zeros(nlayers)
            rt = np.zeros(nlayers)
            if nlayers > 1:
                rt[1:] = rego[y,x][-nlayers+1:]
                rt[0] = bothickness
                mt = melt[y,x][-nlayers:]
                mt[0] = mt[0] * (rt[0]/rego[y,x][-nlayers])
                mf[y,x] = np.sum(mt) / (depth * pix * pix)
            else:
                mf[y,x] = ((melt[y,x][-1] * (depth/rego[y,x][-1])) / (depth * pix * pix))
            
    f1 = plt.figure(1)
    a1 = f1.add_subplot(111)
    if not log:
        md = a1.imshow(mf,origin='lower')#,vmin=0.0,vmax=1.0)
    else:
        md = a1.imshow(mf,origin='lower', norm=colors.LogNorm())#vmin=1e-3))#,vmax=1.0))
    cbar = f1.colorbar(md, label='Total melt fraction')
    if apollo:
        a1.scatter(apollo_x, apollo_y, color='red', label='Apollo 14-17 Landing Sites')
        plt.legend()
    plt.title('Total melt fraction to a depth of ' + str(depth) + ' m')
    plt.show()

    return

if __name__ == '__main__':
    path = '/Users/owner/Documents/git/CTEM/examples/global-lunar-bombardment/'
    #path = '/Users/owner/Desktop/debug/mftest/'
    gridsize = 500
    pix = 12.32e3
    pixkm = pix / 1000

    apollo_x, apollo_y = find_apollo_coords(gridsize,pix)
    surfcoords = np.zeros(len(apollo_y))
    for i in range(len(apollo_y)):
        surfcoords[i] = gridsize - apollo_y[i]

    rego = ctem.util.read_linked_list_binary(path + '/surface_rego.dat', gridsize)
    melt = ctem.util.read_linked_list_binary(path + '/surface_melt.dat', gridsize)
    stack = ctem.util.read_unformatted_binary(path + '/surface_stacknum.dat', gridsize, kind='I4B')
    meltdist = ctem.util.read_linked_list_binary(path + '/surface_meltdist.dat', gridsize, kind='SP')

    surfimg = path + '/surf/surf000001.png'
    img = mpimg.imread(surfimg)

    #nlist = 117
    # craters = ["South Pole-Aitken", "Fecunditatis", "Australe North", "TOPO-13", "Medii", "Bartels-Voskresenskiy",
    #            "Serenitatis North", "Lamont", "Vaporum", "Fowler-Charlier", "Mutus-Vlacq", "Asperitatis", "Aestuum",
    #            "Copernicus-H", "Balmer-Kapetyn", "Cruger-Sirsalis", "Dirichlet-Jackson", "Coloumb-Sarton",
    #            "Schiller-Zucchius", "Smythii", "Amundsen-Ganswindt", "Nubium", "Rupes Recta", "Deslandres", "Lorentz",
    #            "Pozcobutt", "Fitzgerald-Jackson", "Poincare", "Wegener-Winlock", "Orientale Southwest", "TOPO-22",
    #            "Galois", "Schickard", "Fermi", "Ingenii", "Gagarin", "Keeler West", "Birkhoff", "Harkhebi", "Campbell",
    #            "Serenitatis", "Freundlich-Sharanov", "Landau", "Von Karman M", "Leibnitz", "Apollo", "Milne",
    #            "Grimaldi", "Pasteur", "Nectaris", "Mendel-Rydberg", "Moscoviense North", "Moscoviense", "Korolev",
    #            "Mendeleev", "Planck", "Hertzsprung", "Humorum", "Crisisum", "Humboldtianum", "Sikorsky-Rittenhaus",
    #            "Bel'Kovich", "Bailly", "d'Alembert", "Oppenheimer", "Clavius", "Crisium East", "SChwarzschild",
    #            "La Condamine", "Maupertuis", "Imbrium", "Schrodinger", "Ukert", "Campanus", "Protagoras", "Calippus",
    #            "Babbage A", "Prinz", "Lambert R", "Egede", "Fontenelle", "T. Mayer", "Stadius", "Cassini", "Letronne",
    #            "Cruger", "Marius", "Arzachel", "T. Mayer W", "Wallace", "Kepler D", "Orientale", "Lansberg", "Humboldt",
    #            "Kies", "Timaeus", "Bonpland D", "Krieger", "Lansberg C", "Briggs", "Archimedes", "Thebit", "Hansteen",
    #            "Davy", "Markov", "Iridum", "Seleucus", "Cardanus", "Krafft", "Dalton", "St. George", "Bianchinni",
    #            "Herotodus", "Plato", "Damoiseau", "Mairan", "Local"]

    nlist = 8
    craters = ['South Pole-Aitken','Serenitatis','Nectaris','Crisium','Imbrium','Orientale','Iridum','Local']

    aggregate(apollo_x[2],apollo_y[2],5,0,rego,melt,meltdist)
