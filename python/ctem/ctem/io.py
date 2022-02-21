import numpy as np
import os
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LightSource
import matplotlib.cm as cm


#Set pixel scaling common for image writing, at 1 pixel/ array element
dpi = 72.0

def real2float(realstr):
    """
    Converts a Fortran-generated ASCII string of a real value into a np float type. Handles cases where double precision
    numbers in exponential notation use 'd' or 'D' instead of 'e' or 'E'

    Parameters
    ----------
    realstr : string
        Fortran-generated ASCII string of a real value.

    Returns
    -------
     : float
        The converted floating point value of the input string
    """
    return float(realstr.replace('d', 'E').replace('D', 'E'))

def read_param(parameters):
    # Read and parse ctem.in file
    inputfile = os.path.join(parameters['workingdir'], parameters['ctemfile'])

    # Read ctem.in file
    print('Reading input file ' + parameters['ctemfile'])
    with open(inputfile, 'r') as fp:
        lines = fp.readlines()

    # Process file text
    for line in lines:
        fields = line.split()
        if len(fields) > 0:
            if ('pix' == fields[0].lower()): parameters['pix'] = real2float(fields[1])
            if ('gridsize' == fields[0].lower()): parameters['gridsize'] = int(fields[1])
            if ('seed' == fields[0].lower()): parameters['seed'] = int(fields[1])
            if ('sfdfile' == fields[0].lower()): parameters['sfdfile'] = fields[1]
            if ('impfile' == fields[0].lower()): parameters['impfile'] = fields[1]
            if ('maxcrat' == fields[0].lower()): parameters['maxcrat'] = real2float(fields[1])
            if ('sfdcompare' == fields[0].lower()): parameters['sfdcompare'] = fields[1]
            if ('interval' == fields[0].lower()): parameters['interval'] = real2float(fields[1])
            if ('numintervals' == fields[0].lower()): parameters['numintervals'] = int(fields[1])
            if ('popupconsole' == fields[0].lower()): parameters['popupconsole'] = fields[1]
            if ('saveshaded' == fields[0].lower()): parameters['saveshaded'] = fields[1]
            if ('saverego' == fields[0].lower()): parameters['saverego'] = fields[1]
            if ('savepres' == fields[0].lower()): parameters['savepres'] = fields[1]
            if ('savetruelist' == fields[0].lower()): parameters['savetruelist'] = fields[1]
            if ('runtype' == fields[0].lower()): parameters['runtype'] = fields[1]
            if ('restart' == fields[0].lower()): parameters['restart'] = fields[1]
            if ('shadedminh' == fields[0].lower()):
                parameters['shadedminh'] = real2float(fields[1])
                parameters['shadedminhdefault'] = 0
            if ('shadedmaxh' == fields[0].lower()):
                parameters['shadedmaxh'] = real2float(fields[1])
                parameters['shadedmaxhdefault'] = 0

    # Test values for further processing
    if (parameters['interval'] <= 0.0):
        print('Invalid value for or missing variable INTERVAL in ' + inputfile)
    if (parameters['numintervals'] <= 0):
        print('Invalid value for or missing variable NUMINTERVALS in ' + inputfile)
    if (parameters['pix'] <= 0.0):
        print('Invalid value for or missing variable PIX in ' + inputfile)
    if (parameters['gridsize'] <= 0):
        print('Invalid value for or missing variable GRIDSIZE in ' + inputfile)
    if (parameters['seed'] == 0):
        print('Invalid value for or missing variable SEED in ' + inputfile)
    if (parameters['sfdfile'] is None):
        print('Invalid value for or missing variable SFDFILE in ' + inputfile)
    if (parameters['impfile'] is None):
        print('Invalid value for or missing variable IMPFILE in ' + inputfile)
    if (parameters['popupconsole'] is None):
        print('Invalid value for or missing variable POPUPCONSOLE in ' + inputfile)
    if (parameters['saveshaded'] is None):
        print('Invalid value for or missing variable SAVESHADED in ' + inputfile)
    if (parameters['saverego'] is None):
        print('Invalid value for or missing variable SAVEREGO in ' + inputfile)
    if (parameters['savepres'] is None):
        print('Invalid value for or missing variable SAVEPRES in ' + inputfile)
    if (parameters['savetruelist'] is None):
        print('Invalid value for or missing variable SAVETRUELIST in ' + inputfile)
    if (parameters['runtype'] is None):
        print('Invalid value for or missing variable RUNTYPE in ' + inputfile)
    if (parameters['restart'] is None):
        print('Invalid value for or missing variable RESTART in ' + inputfile)
    
    return parameters


def read_formatted_ascii(filename, skip_lines):
    # Generalized ascii text reader
    # For use with sfdcompare, production, odist, tdist, pdist data files
    data = np.genfromtxt(filename, skip_header=skip_lines)
    return data


def read_unformatted_binary(filename, gridsize):
    # Read unformatted binary files created by Fortran
    # For use with surface ejecta and surface dem data files
    dt = np.float
    data = np.fromfile(filename, dtype=dt)
    data.shape = (gridsize, gridsize)
    
    return data


def read_ctemdat(parameters, seedarr):
    # Read and parse ctem.dat file
    datfile = parameters['workingdir'] + 'ctem.dat'

    # Read ctem.dat file
    print('Reading input file ' + parameters['datfile'])
    fp = open(datfile, 'r')
    lines = fp.readlines()
    fp.close()

    # Parse file lines and update parameter fields
    fields = lines[0].split()
    if len(fields) > 0:
        parameters['totalimpacts'] = real2float(fields[0])
        parameters['ncount'] = int(fields[1])
        parameters['curyear'] = real2float(fields[2])
        parameters['restart'] = fields[3]
        parameters['fracdone'] = real2float(fields[4])
        parameters['masstot'] = real2float(fields[5])

    # Parse remainder of file to build seed array
    nlines = len(lines)
    index = 1
    while (index < nlines):
        fields = lines[index].split()
        seedarr[index - 1] = real2float(fields[0])
        index += 1
    
    parameters['seedn'] = index - 1
    
    return


def read_impact_mass(filename):
    # Read impact mass file
    
    fp = open(filename, 'r')
    line = fp.readlines()
    fp.close()
    
    fields = line[0].split()
    if (len(fields) > 0):
        mass = real2float(fields[0])
    else:
        mass = 0
    
    return mass


# Write production function to file production.dat
# This file format does not exactly match that generated from IDL. Does it work?
def write_production(parameters, production):
    filename = parameters['workingdir'] + parameters['sfdfile']
    np.savetxt(filename, production, fmt='%1.8e', delimiter='   ')
    
    return


def create_dir_structure(parameters):
    # Create directories for various output files if they do not already exist
    directories = ['dist', 'misc', 'rego', 'rplot', 'surf', 'shaded']
    
    for directory in directories:
        dir_test = parameters['workingdir'] + directory
        if not os.path.isdir(dir_test):
            os.makedirs(dir_test)
    
    return


def write_ctemdat(parameters, seedarr):
    # Write various parameters and random number seeds into ctem.dat file
    filename = parameters['workingdir'] + parameters['datfile']
    fp = open(filename, 'w')
    
    template = "%(totalimpacts)17d %(ncount)12d %(curyear)19.12E %(restart)s %(fracdone)9.6f %(masstot)19.12E\n"
    fp.write(template % parameters)

    # Write random number seeds to the file
    for index in range(parameters['seedn']):
        fp.write("%12d\n" % seedarr[index])
    
    fp.close()
    
    return


def copy_dists(parameters):
    # Save copies of distribution files
    
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


# Possible references
# http://nbviewer.jupyter.org/github/ThomasLecocq/geophysique.be/blob/master/2014-02-25%20Shaded%20Relief%20Map%20in%20Python.ipynb   

def image_dem(parameters, DEM):
    dpi = 300.0  # 72.0
    pix = parameters['pix']
    gridsize = parameters['gridsize']
    ve = 1.0
    azimuth = 300.0  # parameters['azimuth']
    solar_angle = 20.0  # parameters['solar_angle']
    
    ls = LightSource(azdeg=azimuth, altdeg=solar_angle)
    dem_img = ls.hillshade(DEM, vert_exag=ve, dx=pix, dy=pix)
    
    # Generate image to put into an array
    height = gridsize / dpi
    width = gridsize / dpi
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    ax = plt.axes([0, 0, 1, 1])
    ax.imshow(dem_img, interpolation="nearest", cmap='gray', vmin=0.0, vmax=1.0)
    plt.axis('off')
    # Save image to file
    filename = parameters['workingdir'] + 'surf' + os.sep + "surf%06d.png" % parameters['ncount']
    plt.savefig(filename, dpi=dpi, bbox_inches=0)
    
    return


def image_shaded_relief(parameters, DEM):
    dpi = 300.0  # 72.0
    pix = parameters['pix']
    gridsize = parameters['gridsize']
    ve = 1.0
    mode = 'overlay'
    azimuth = 300.0  # parameters['azimuth']
    solar_angle = 20.0  # parameters['solar_angle']
    
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
        shadedminh = np.amin(DEM)
    else:
        shadedminh = parameters['shadedminh']
    if (parameters['shadedmaxhdefault'] == 1):
        shadedmaxh = np.amax(DEM)
    else:
        shadedmaxh = parameters['shadedmaxh']
    
    dem_img = ls.shade(DEM, cmap=cmap, blend_mode=mode, fraction=1.0,
                       vert_exag=ve, dx=pix, dy=pix,
                       vmin=shadedminh, vmax=shadedmaxh)

    # Generate image to put into an array
    height = gridsize / dpi
    width = gridsize / dpi
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    ax = plt.axes([0, 0, 1, 1])
    ax.imshow(dem_img, interpolation="nearest", vmin=0.0, vmax=1.0)
    plt.axis('off')
    # Save image to file
    filename = parameters['workingdir'] + 'shaded' + os.sep + "shaded%06d.png" % parameters['ncount']
    plt.savefig(filename, dpi=dpi, bbox_inches=0)
    return parameters


def image_regolith(parameters, regolith):
    # Create scaled regolith image
    minref = parameters['pix'] * 1.0e-4
    maxreg = np.amax(regolith)
    minreg = np.amin(regolith)
    if (minreg < minref): minreg = minref
    if (maxreg < minref): maxreg = (minref + 1.0e3)
    regolith_scaled = np.copy(regolith)
    np.place(regolith_scaled, regolith_scaled < minref, minref)
    regolith_scaled = 254.0 * (
            (np.log(regolith_scaled) - np.log(minreg)) / (np.log(maxreg) - np.log(minreg)))
    
    # Save image to file
    filename = parameters['workingdir'] + 'rego' + os.sep + "rego%06d.png" % parameters['ncount']
    height = parameters['gridsize'] / dpi
    width = height
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    fig.figimage(regolith_scaled, cmap=cm.nipy_spectral, origin='lower')
    plt.savefig(filename)
    
    return


def create_rplot(parameters, odist, pdist, tdist, ph1):
    # Parameters: empirical saturation limit and dfrac
    satlimit = 3.12636
    dfrac = 2 ** (1. / 4) * 1.0e-3

    # Calculate geometric saturation
    minx = (parameters['pix'] / 3.0) * 1.0e-3
    maxx = 3 * parameters['pix'] * parameters['gridsize'] * 1.0e-3
    geomem = np.array([[minx, satlimit / 20.0], [maxx, satlimit / 20.0]])
    geomep = np.array([[minx, satlimit / 10.0], [maxx, satlimit / 10.0]])

    # Create distribution arrays without zeros for plotting on log scale
    idx = np.nonzero(odist[:, 5])
    odistnz = odist[idx]
    odistnz = odistnz[:, [2, 5]]
    
    idx = np.nonzero(pdist[:, 5])
    pdistnz = pdist[idx]
    pdistnz = pdistnz[:, [2, 5]]
    
    idx = np.nonzero(tdist[:, 5])
    tdistnz = tdist[idx]
    tdistnz = tdistnz[:, [2, 5]]

    # Correct pdist
    pdistnz[:, 1] = pdistnz[:, 1] * parameters['curyear'] / parameters['interval']

    # Create sdist bin factors, which contain one crater per bin
    area = (parameters['gridsize'] * parameters['pix'] * 1.0e-3) ** 2.
    plo = 1
    sq2 = 2 ** (1. / 2)
    while (sq2 ** plo > minx):
        plo = plo - 1
    phi = plo + 1
    while (sq2 ** phi < maxx):
        phi = phi + 1
    n = phi - plo + 1
    sdist = np.zeros([n, 2])
    p = plo
    for index in range(n):
        sdist[index, 0] = sq2 ** p
        sdist[index, 1] = sq2 ** (2.0 * p + 1.5) / (area * (sq2 - 1))
        p = p + 1

    # Create time label
    tlabel = "%5.4e" % parameters['curyear']
    tlabel = tlabel.split('e')
    texp = str(int(tlabel[1]))
    timelabel = 'Time = ' + r'${}$ x 10$^{}$'.format(tlabel[0], texp) + ' yrs'

    # Save image to file
    filename = parameters['workingdir'] + 'rplot' + os.sep + "rplot%06d.png" % parameters['ncount']
    height = parameters['gridsize'] / dpi
    width = height
    fig = plt.figure(figsize=(width, height), dpi=dpi)

    # Alter background color to be black, and change axis colors accordingly
    plt.style.use('dark_background')
    plt.rcParams['axes.prop_cycle']

    # Plot data    
    plt.plot(odistnz[:, 0] * 1.0e-3, odistnz[:, 1], linewidth=3.0, color='blue')
    plt.plot(pdistnz[:, 0] * 1.0e-3, pdistnz[:, 1], linewidth=2.0, linestyle='dashdot', color='white')
    plt.plot(tdistnz[:, 0] * 1.0e-3, tdistnz[:, 1], linewidth=2.0, color='red')
    
    plt.plot(geomem[:, 0], geomem[:, 1], linewidth=2.0, linestyle=':', color='yellow')
    plt.plot(geomep[:, 0], geomep[:, 1], linewidth=2.0, linestyle=':', color='yellow')
    plt.plot(sdist[:, 0], sdist[:, 1], linewidth=2.0, linestyle=':', color='yellow')
    
    plt.plot(ph1[:, 0] * dfrac, ph1[:, 1], 'wo')

    # Create plot labels    
    plt.title('Crater Distribution R-Plot', fontsize=22)
    plt.xlim([minx, maxx])
    plt.xscale('log')
    plt.xlabel('Crater Diameter (km)', fontsize=18)
    plt.ylim([5.0e-4, 5.0])
    plt.yscale('log')
    plt.ylabel('R Value', fontsize=18)
    plt.text(1.0e-2, 1.0, timelabel, fontsize=18)
    
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=12)
    
    plt.savefig(filename)
    
    return
