#!/usr/local/bin/python
#
#Cratered Terrain Evolution Model file and graphical writing utilities
#
#Original IDL design: Jim Richardson, Arecibo Observatory
#Revised IDL design: David Minton, Purdue University
#Re-engineered design and Python implementation: Matthew Route, Purdue University
#August 2016
#

import matplotlib
from matplotlib import pyplot
import numpy
import os
import shutil
import scipy
from scipy import signal
from matplotlib.colors import LightSource
import matplotlib.cm as cm

#Set pixel scaling common for image writing, at 1 pixel/ array element
dpi = 72.0

#Write production function to file production.dat
#This file format does not exactly match that generated from IDL. Does it work?
def write_production(parameters, production):
    filename = parameters['workingdir'] + parameters['sfdfile']
    numpy.savetxt(filename, production, fmt='%1.8e', delimiter='   ')
    
    return

def create_dir_structure(parameters):
    #Create directories for various output files if they do not already exist
    directories=['dist','misc','rego','rplot','surf','shaded']
    
    for directory in directories:
        dir_test = parameters['workingdir'] + directory
        if not os.path.isdir(dir_test):
            os.makedirs(dir_test)
    
    return

def write_ctemdat(parameters, seedarr):
    #Write various parameters and random number seeds into ctem.dat file
    filename = parameters['workingdir'] + parameters['datfile']
    fp = open(filename,'w')
    
    template = "%(totalimpacts)17d %(ncount)12d %(curyear)19.12E %(restart)s %(fracdone)9.6f %(masstot)19.12E\n"
    fp.write(template % parameters)
    
    #Write random number seeds to the file
    for index in range(parameters['seedn']):
       fp.write("%12d\n" % seedarr[index])
       
    fp.close()   
    
    return

def copy_dists(parameters):
    #Save copies of distribution files
        
    orig_list = ['odistribution', 'ocumulative', 'pdistribution', 'tdistribution']
    dest_list = ['odist', 'ocum', 'pdist', 'tdist']       
        
    for index in range(len(orig_list)):  
        forig = parameters['workingdir'] + orig_list[index] + '.dat'
        fdest = parameters['workingdir'] + 'dist' + os.sep + dest_list[index] + "_%06d.dat" % parameters['ncount']
        shutil.copy2(forig, fdest)
    
    forig = parameters['workingdir'] + 'impactmass.dat'
    fdest = parameters['workingdir'] + 'misc' + os.sep + "mass_%06d.dat" % parameters['ncount']
    shutil.copy2(forig, fdest)   
    
    if (parameters['savetruelist'].upper() == 'T'):
        forig = parameters['workingdir'] + 'tcumulative.dat'
        fdest = parameters['workingdir'] + 'dist' + os.sep + "tcum_%06d.dat" % parameters['ncount']
        shutil.copy2(forig, fdest)
    
    return

#Possible references
#http://nbviewer.jupyter.org/github/ThomasLecocq/geophysique.be/blob/master/2014-02-25%20Shaded%20Relief%20Map%20in%20Python.ipynb   
   
def image_dem(parameters, DEM):
    dpi = 300.0  # 72.0
    pix = parameters['pix']
    gridsize = parameters['gridsize']
    ve = 1.0
    azimuth = 300.0 #parameters['azimuth']
    solar_angle = 20.0 #parameters['solar_angle']

    ls = LightSource(azdeg=azimuth, altdeg=solar_angle)
    dem_img = ls.hillshade(DEM, vert_exag=ve, dx=pix, dy=pix)

    # Generate image to put into an array
    height = gridsize / dpi
    width = gridsize / dpi
    fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=dpi)
    ax = matplotlib.pyplot.axes([0, 0, 1, 1])
    ax.imshow(dem_img, interpolation="nearest", cmap='gray',vmin=0.0,vmax=1.0)
    matplotlib.pyplot.axis('off')
    # Save image to file
    filename = parameters['workingdir'] + 'surf' + os.sep + "surf%06d.png" % parameters['ncount']
    matplotlib.pyplot.savefig(filename,dpi=dpi,bbox_inches=0)

    return

def image_shaded_relief(parameters, DEM):
    dpi = 300.0 #72.0
    pix = parameters['pix']
    gridsize = parameters['gridsize']
    ve = 1.0
    mode = 'overlay'
    azimuth = 300.0 #parameters['azimuth']
    solar_angle = 20.0  #parameters['solar_angle']

    ls = LightSource(azdeg=azimuth, altdeg=solar_angle)
    cmap = cm.cividis

    # If min and max appear to be reversed, then fix them
    if (parameters['shadedminh'] > parameters['shadedmaxh']):
        temp = parameters['shadedminh']
        parameters['shadedminh'] = parameters['shadedmaxh']
        parameters['shadedmaxh'] = temp
    else:
        parameters['shadedminh'] = parameters['shadedminh']
        parameters['shadedmaxh'] = parameters['shadedmaxh']

    # If no shadedmin/max parameters are read in from ctem.dat, determine the values from the data
    if (parameters['shadedminhdefault'] == 1):
        shadedminh = numpy.amin(DEM)
    else:
        shadedminh = parameters['shadedminh']
    if (parameters['shadedmaxhdefault'] == 1):
        shadedmaxh = numpy.amax(DEM)
    else:
        shadedmaxh = parameters['shadedmaxh']

    dem_img = ls.shade(DEM, cmap=cmap,blend_mode=mode, fraction=1.0,
                       vert_exag=ve, dx=pix, dy=pix,
                       vmin=shadedminh, vmax=shadedmaxh)

    #Generate image to put into an array
    height = gridsize / dpi
    width = gridsize / dpi
    fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=dpi)
    ax = matplotlib.pyplot.axes([0, 0, 1, 1])
    ax.imshow(dem_img,interpolation="nearest",vmin=0.0,vmax=1.0)
    matplotlib.pyplot.axis('off')
    # Save image to file
    filename = parameters['workingdir'] + 'shaded' + os.sep + "shaded%06d.png" % parameters['ncount']
    matplotlib.pyplot.savefig(filename,dpi=dpi,bbox_inches=0)
    return


def image_regolith(parameters, regolith):
    # Create scaled regolith image
    minref = parameters['pix'] * 1.0e-4
    maxreg = numpy.amax(regolith)
    minreg = numpy.amin(regolith)
    if (minreg < minref): minreg = minref
    if (maxreg < minref): maxreg = (minref + 1.0e3)
    regolith_scaled = numpy.copy(regolith)
    numpy.place(regolith_scaled, regolith_scaled < minref, minref)
    regolith_scaled = 254.0 * (
                (numpy.log(regolith_scaled) - numpy.log(minreg)) / (numpy.log(maxreg) - numpy.log(minreg)))

    # Save image to file
    filename = parameters['workingdir'] + 'rego' + os.sep + "rego%06d.png" % parameters['ncount']
    height = parameters['gridsize'] / dpi
    width = height
    fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=dpi)
    fig.figimage(regolith_scaled, cmap=matplotlib.cm.nipy_spectral, origin='lower')
    matplotlib.pyplot.savefig(filename)

    return

def create_rplot(parameters,odist,pdist,tdist,ph1):
    #Parameters: empirical saturation limit and dfrac
    satlimit = 3.12636
    dfrac = 2**(1./4) * 1.0e-3
    
    #Calculate geometric saturation
    minx = (parameters['pix'] / 3.0) * 1.0e-3 
    maxx = 3 * parameters['pix'] * parameters['gridsize'] * 1.0e-3     
    geomem = numpy.array([[minx, satlimit / 20.0], [maxx, satlimit / 20.0]])
    geomep = numpy.array([[minx, satlimit / 10.0], [maxx, satlimit / 10.0]])
    
    #Create distribution arrays without zeros for plotting on log scale
    idx = numpy.nonzero(odist[:,5])
    odistnz = odist[idx]
    odistnz = odistnz[:,[2,5]]
    
    idx = numpy.nonzero(pdist[:,5])
    pdistnz = pdist[idx]
    pdistnz = pdistnz[:,[2,5]]
    
    idx = numpy.nonzero(tdist[:,5])
    tdistnz = tdist[idx]
    tdistnz = tdistnz[:,[2,5]]
    
    #Correct pdist
    pdistnz[:,1] = pdistnz[:,1] * parameters['curyear'] / parameters['interval']    
    
    #Create sdist bin factors, which contain one crater per bin
    area = (parameters['gridsize'] * parameters['pix'] * 1.0e-3)**2.
    plo = 1
    sq2 = 2**(1./2)
    while (sq2**plo > minx):
        plo = plo - 1
    phi = plo + 1
    while (sq2**phi < maxx):
        phi = phi + 1
    n = phi - plo + 1
    sdist = numpy.zeros([n , 2])
    p = plo
    for index in range(n):
        sdist[index, 0] = sq2**p
        sdist[index, 1] = sq2**(2.0*p + 1.5) / (area * (sq2 - 1))
        p = p + 1
    
    #Create time label
    tlabel = "%5.4e" % parameters['curyear']
    tlabel = tlabel.split('e')
    texp = str(int(tlabel[1]))
    timelabel = 'Time = '+ r'${}$ x 10$^{}$'.format(tlabel[0], texp) + ' yrs'    
    
    #Save image to file
    filename = parameters['workingdir'] + 'rplot' + os.sep + "rplot%06d.png" % parameters['ncount']
    height = parameters['gridsize'] / dpi
    width = height
    fig = matplotlib.pyplot.figure(figsize = (width, height), dpi = dpi)
    
    #Alter background color to be black, and change axis colors accordingly
    matplotlib.pyplot.style.use('dark_background')
    matplotlib.pyplot.rcParams['axes.prop_cycle']
    
    #Plot data    
    matplotlib.pyplot.plot(odistnz[:,0]*1.0e-3, odistnz[:,1], linewidth=3.0, color = 'blue')
    matplotlib.pyplot.plot(pdistnz[:,0]*1.0e-3, pdistnz[:,1], linewidth=2.0, linestyle='dashdot', color = 'white')
    matplotlib.pyplot.plot(tdistnz[:,0]*1.0e-3, tdistnz[:,1], linewidth=2.0, color = 'red')
    
    matplotlib.pyplot.plot(geomem[:,0], geomem[:,1], linewidth=2.0, linestyle =':', color = 'yellow')
    matplotlib.pyplot.plot(geomep[:,0], geomep[:,1], linewidth=2.0, linestyle =':', color = 'yellow')
    matplotlib.pyplot.plot(sdist[:,0], sdist[:,1], linewidth=2.0, linestyle =':', color = 'yellow')
    
    matplotlib.pyplot.plot(ph1[:,0] * dfrac, ph1[:,1], 'wo')    
    
    #Create plot labels    
    matplotlib.pyplot.title('Crater Distribution R-Plot',fontsize=22)
    matplotlib.pyplot.xlim([minx, maxx])
    matplotlib.pyplot.xscale('log')
    matplotlib.pyplot.xlabel('Crater Diameter (km)',fontsize=18)
    matplotlib.pyplot.ylim([5.0e-4, 5.0])    
    matplotlib.pyplot.yscale('log')
    matplotlib.pyplot.ylabel('R Value', fontsize=18)
    matplotlib.pyplot.text(1.0e-2, 1.0, timelabel, fontsize=18)
    
    matplotlib.pyplot.tick_params(axis='both', which='major', labelsize=14)
    matplotlib.pyplot.tick_params(axis='both', which='minor', labelsize=12)
    
    matplotlib.pyplot.savefig(filename)   
   
    return