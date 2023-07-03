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

def find_all_ages(x,y,depth,rego,age):
    nlayers, bothickness = traverse(y,x,depth,rego)
    agearray = np.zeros((nlayers,nage))
    rt = np.zeros(nlayers)
    rt[0] = bothickness
    if nlayers > 1:
        rt[1:] = rego[y,x][-nlayers+1:]
    for i in range(1,nlayers+1):
        if i == 1:
            agearray[i-1] = age[y,x][-nage:]
        else:
            agearray[i-1] = age[y,x][-(i*nage):-(nage*(i-1))]
    agearray[nlayers-1] = agearray[nlayers-1] * (rt[0]/rego[y,x][-nlayers])
    
    agedist = np.zeros(nage)
    for j in range(nage):
        agedist[j] = np.sum(agearray[:,j])
        
    
    trueages = []
    xax = np.arange(nage)
    for a in xax:
        trueage = ctem.craterproduction.T_from_scale(((1+a)*(interval/nage)),'NPF_Moon')
        trueages.append(trueage[0])
    order = trueages[::-1]
    
    totalmelt, idc = find_melt_at_pixel(x,y,depth,rego,melt,meltdist)
    qmcarray = np.zeros(nage)
    for i in range(nlist-1):
        age_number = qmcages[i]
        for j in range(nage):
            if age_number > order[j]:
                qmcarray[j] = qmcarray[j] + totalmelt[i]
                break
            elif j == nage-1:
                qmcarray[j] = qmcarray[j] + totalmelt[i]
            else:
                continue
                
    return order, agedist, qmcarray

def aggregate_ages(x,y,depth,n,rego,age,log=False,apollo=False):
    if n==0:
        order,agedist,qmcarray = find_all_ages(x,y,depth,rego,age)
        k = 1
    elif n<0:
        print("number of pixels must be >=0")
    else:
        agedist = []
        qmcarray = []
        k = 0
        for i in range(-n,n+1):
            for j in range(-n,n+1):
                agedist.append(np.zeros(nage))
                qmcarray.append(np.zeros(nage))
                o,localage,qmcage = find_all_ages(x+i,y+j,depth,rego,age)
                agedist[k] = localage
                qmcarray[k] = qmcage
                k += 1
        order = o
        
    agedist = np.asarray(agedist)
    qmcarray = np.asarray(qmcarray)

    totl = np.zeros(nage)
    totq = np.zeros(nage)

    for l in range(nage):
        if n !=0:
            totl[l] = np.sum(agedist[:,l])
            totq[l] = np.sum(qmcarray[:,l])
        else:
            totl[l] = agedist[l]
            totq[l] = qmcarray[l]
            
    totl = totl / (k * pix * pix)
    totq = totq / (k * pix * pix)
    
    mf = (np.sum(totl) + np.sum (totq)) / depth
                
    
    #plot
    
    f1, (a1, a2) = plt.subplots(1,2)
    f1.suptitle('Pixel: [' + '{:3d}'.format(x) + ',' + '{:3d}'.format(y) + '] (Aggregated to ' + str(n) + ' pixels) (Total melt fraction = ' + '{0:.4f}'.format(mf) + ') (depth = ' + '{0:.2f}'.format(depth) + ' m)')
    a1.bar(order,totl,width=((order[0]-order[-1])/nage),log=log,color='#1f78b4',label='Local Melt')
    a1.bar(order,totq,width=((order[0]-order[-1])/nage),bottom=totl,log=log,color='#d95f02',label='Quasi-MC Melt')
    a1.set_xlabel('Age (Ga)')
    a1.set_ylabel('Melt Depth (m)')
    a2.imshow(img)
    if n == 0:
        a2.scatter(x, gridsize-y, color='blue', label='Selected Region')
    else:
        rect = patches.Rectangle((x,gridsize-y),n,n,linewidth=1,edgecolor='b',facecolor='b',label='Selected Region')
        a2.add_patch(rect)
    a2.xaxis.set_ticklabels([])
    a2.yaxis.set_ticklabels([])
    if apollo:
        a2.scatter(apollo_x, surfcoords, color='red', label='Apollo 14-17 Landing Sites')
    a2.tick_params(left = False)
    a2.tick_params(bottom = False)
    a1.legend()
    a2.legend()
    plt.rcParams['figure.figsize'] = [13, 5]
    plt.show()
    return mf * depth

if __name__ == '__main__':
    path = '/Users/owner/Documents/git/CTEM/examples/global-lunar-bombardment/'
    #path = '/Users/owner/Desktop/debug/mftest/'
    #gridsize = 500
    #pix = 12.32e3
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
    age = ctem.util.read_linked_list_binary(path + '/surface_age.dat', gridsize, kind='SP')


    surfimg = path + '/surf/surf000001.png'
    img = mpimg.imread(surfimg)

    interval = 628.0
    nage = 60

    nlist = 117
    craters = ["South Pole-Aitken", "Fecunditatis", "Australe North", "TOPO-13", "Medii", "Bartels-Voskresenskiy",
               "Serenitatis North", "Lamont", "Vaporum", "Fowler-Charlier", "Mutus-Vlacq", "Asperitatis", "Aestuum",
               "Copernicus-H", "Balmer-Kapetyn", "Cruger-Sirsalis", "Dirichlet-Jackson", "Coloumb-Sarton",
               "Schiller-Zucchius", "Smythii", "Amundsen-Ganswindt", "Nubium", "Rupes Recta", "Deslandres", "Lorentz",
               "Pozcobutt", "Fitzgerald-Jackson", "Poincare", "Wegener-Winlock", "Orientale Southwest", "TOPO-22",
               "Galois", "Schickard", "Fermi", "Ingenii", "Gagarin", "Keeler West", "Birkhoff", "Harkhebi", "Campbell",
               "Serenitatis", "Freundlich-Sharanov", "Landau", "Von Karman M", "Leibnitz", "Apollo", "Milne",
               "Grimaldi", "Pasteur", "Nectaris", "Mendel-Rydberg", "Moscoviense North", "Moscoviense", "Korolev",
               "Mendeleev", "Planck", "Hertzsprung", "Humorum", "Crisisum", "Humboldtianum", "Sikorsky-Rittenhaus",
               "Bel'Kovich", "Bailly", "d'Alembert", "Oppenheimer", "Clavius", "Crisium East", "SChwarzschild",
               "La Condamine", "Maupertuis", "Imbrium", "Schrodinger", "Ukert", "Campanus", "Protagoras", "Calippus",
               "Babbage A", "Prinz", "Lambert R", "Egede", "Fontenelle", "T. Mayer", "Stadius", "Cassini", "Letronne",
               "Cruger", "Marius", "Arzachel", "T. Mayer W", "Wallace", "Kepler D", "Orientale", "Lansberg", "Humboldt",
               "Kies", "Timaeus", "Bonpland D", "Krieger", "Lansberg C", "Briggs", "Archimedes", "Thebit", "Hansteen",
               "Davy", "Markov", "Iridum", "Seleucus", "Cardanus", "Krafft", "Dalton", "St. George", "Bianchinni",
               "Herotodus", "Plato", "Damoiseau", "Mairan", "Local"]
    
    qmcages = [4.310, 4.307, 4.304, 4.301, 4.298, 4.295, 4.292, 4.289, 4.286, 4.283, 4.280, 4.277, 4.274, 4.270, 4.265, 4.260, 4.252, 4.249, 4.246, 4.242, 4.239, 4.236, 4.233, 4.231, 4.229, 4.226, 4.223, 4.220, 4.217, 4.214, 4.210, 4.208, 4.206, 4.204, 4.202, 4.200, 4.198, 4.196, 4.194, 4.191, 4.187, 4.183, 4.177, 4.173, 4.169, 4.164, 4.162, 4.160, 4.158, 4.156, 4.152, 4.112, 4.110, 4.105, 4.100, 4.098, 4.095, 4.090, 4.070, 4.048, 3.980, 3.950, 3.940, 3.930, 3.920, 3.910, 3.900, 3.890, 3.880, 3.875, 3.870, 3.860, 3.855, 3.850, 3.846, 3.843, 3.840, 3.838, 3.836, 3.834, 3.832, 3.830, 3.828, 3.826, 3.824, 3.822, 3.820, 3.818, 3.816, 3.814, 3.812, 3.810, 3.800, 3.750, 3.700, 3.650, 3.600, 3.500, 3.500, 3.490, 3.480, 3.470, 3.460, 3.450, 3.420, 3.400, 3.335, 3.333, 3.330, 3.320, 3.310, 3.300, 3.200, 3.100, 3.050, 3.000]

    # nlist = 8
    # craters = ['South Pole-Aitken','Serenitatis','Nectaris','Crisium','Imbrium','Orientale','Iridum','Local']

    #aggregate(apollo_x[2],apollo_y[2],5,0,rego,melt,meltdist)
    aggregate_ages(apollo_x[0],apollo_y[0],5,1,rego,age,log=False,apollo=True)