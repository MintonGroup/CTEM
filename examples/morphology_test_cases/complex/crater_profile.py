import numpy as np
import pandas as pd
import os
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LightSource

def hillshade2(filename,gridsize,pix,array):
    dpi = 300.0 #72.0
    ls = LightSource(azdeg=315, altdeg=45)
    shadimg = ls.shade(array/pix,cmap=cm.gray)

    height = gridsize / dpi
    width = height
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    fig.figimage(shadimg, origin='lower')
    plt.savefig(filename)
    print(f'Image saved to {filename}')
    return

def image_dem(filename,gridsize,pix,CTEM_dem):
    dpi = 300.0 #72.0
    # Create surface dem map
    solar_angle = 20.0
    dem_map = np.copy(CTEM_dem) - np.roll(CTEM_dem, 1, 0)
    dem_map = (0.5 * np.pi) + np.arctan2(dem_map, pix)
    dem_map = dem_map - np.radians(solar_angle) * (0.5 * np.pi)
    np.place(dem_map, dem_map > (0.5 * np.pi), 0.5 * np.pi)
    dem_map = np.absolute(dem_map)
    dem_map = 254.0 * np.cos(dem_map)

    # Save image to file
    height = gridsize / dpi
    width = height
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    fig.figimage(dem_map, cmap=cm.gray, origin='lower')
    plt.savefig(filename)
    print(f'Image saved to {filename}')

    return


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

crater_name = "Copernicus"
path = f'CSV/{crater_name}.csv'


unit_correction = 0.2727 #Ratio of Earth to Moon radii
crater = pd.read_csv(path) #  * unit_correction * 1e-3
crater['X'] = crater['X'] * unit_correction * 1e-3
crater['Y'] = crater['Y'] * unit_correction * 1e-3
crater['Z'] = crater['Z'] * np.sqrt(unit_correction) * 1e-3


#print('Original input file head')
#print(crater.head())
#
xshift = (crater['X'].max() + crater['X'].min())/2
yshift = (crater['Y'].max() + crater['Y'].min())/2
zshift = -1.0
crater['X'] = crater['X'] - xshift
crater['Y'] = crater['Y'] - yshift
crater['Z'] = crater['Z'] - zshift
crater = crater.sort_values(by=['Y','X'])

#
#print('Centered and sorted values head')
#print(crater.head())
#
gridsize = int(np.sqrt(crater.shape[0]))
#
#print(f'Box size is {gridsize} pixels')
#
xres = 2 * crater['X'].max() / gridsize
yres = 2 * crater['Y'].max() / gridsize
#zres = 0.5 * (xres + yres)
print('X and Y resolutions should be the same')
print(f'X resolution is {xres} km / pix')
print(f'Y resolution is {yres} km / pix')

xdata = crater['X'].values * xres
ydata = crater['Y'].values * xres
#
pix = xres
xdata = crater['X'].values.reshape(gridsize,gridsize)[0,:]
ydata = crater['Y'].values.reshape(gridsize,gridsize)[:,0]
orig_dem = crater['Z'].values.reshape(gridsize,gridsize)

f = interpolate.interp2d(xdata, ydata, orig_dem, kind='cubic')

proffig = plt.figure(1, figsize=(10,2))
axes = proffig.add_subplot(111)
axes.set_xlabel("Distance from center (km)")
axes.set_ylabel("Elevantion (km)")

rho = np.arange(0.0,xdata.max(),step=pix)
phi = np.arange(0,2*np.pi,step=30*np.pi/360.0)
zinterp = np.empty([rho.size,phi.size])

xprof = np.empty_like(zinterp)
yprof = np.empty_like(zinterp)

for jdx,p in enumerate(phi):
    for idx, r in enumerate(rho):
        x,y = pol2cart(r,p)
        xprof[idx,jdx], yprof[idx,jdx] = pol2cart(r, p)
        zinterp[idx,jdx] = f(x,y)
        xprof[idx,jdx] = x
        yprof[idx, jdx] = y
    axes.plot(rho,zinterp[:,jdx],':')


#Compare with CTEM output
gridsize = 2000
dem_file = 'surface_dem.dat'
CTEM_dem = np.fromfile(dem_file , dtype = np.float64)
CTEM_dem.shape = (gridsize,gridsize)
CTEM_dem *= 1e-3
pix = 0.200
gridsize = 2000
CTEM_dem.shape = (gridsize,gridsize)
xdata, ydata = np.mgrid[-gridsize/2:gridsize/2, -gridsize/2:gridsize/2] * pix

xdata = np.arange(-gridsize/2,gridsize/2,1) * pix + pix/2
ydata = np.arange(-gridsize/2,gridsize/2,1) * pix + pix/2

f = interpolate.interp2d(xdata, ydata, CTEM_dem, kind='cubic')

rho = np.arange(0.0,xdata.max(),step=pix)
phi = np.arange(0,2*np.pi,step=30*np.pi/360.0)
zinterp = np.empty([rho.size,phi.size])

xprof = np.empty_like(zinterp)
yprof = np.empty_like(zinterp)

for jdx,p in enumerate(phi):
    for idx, r in enumerate(rho):
        x,y = pol2cart(r,p)
        xprof[idx,jdx], yprof[idx,jdx] = pol2cart(r, p)
        zinterp[idx,jdx] = f(x,y)
        xprof[idx,jdx] = x
        yprof[idx, jdx] = y
    axes.plot(rho,zinterp[:,jdx],'-')


plt.savefig(f'{crater_name}_profile.png',dpi=300,bbox_inches='tight')
os.system(f'open {crater_name}_profile.png')

imgfile = f'{crater_name}.png'
image_dem(imgfile,gridsize,pix,orig_dem)
#os.system(f'open {imgfile}')

#Let's reconstruct the original DEM by interpolating betwen profiles
#xp = xprof.flatten()
#yp = yprof.flatten()
#zp = zinterp.flatten()
#
#grid_x, grid_y = np.mgrid[-gridsize/2:gridsize/2, -gridsize/2:gridsize/2] * pix
#points = (xp,yp)
#
#zrecon = interpolate.griddata(points,zp,(grid_x,grid_y),method='cubic')
#
#imgfile = f'{crater_name}-reconstructed.png'
#image_dem(imgfile,gridsize,pix,zrecon)
#os.system(f'open {imgfile}')

