import numpy as np
import matplotlib.pyplot as plt
import ctem
import os

#path = '/Users/owner/Desktop/yh/ctem_20220808'
#path = '/Users/owner/Desktop/raytest/off'
#path = '/Users/owner/Desktop/melt_distribution/'
path = os.getcwd()
gridsize = 2000
pix = 2000
pixkm = pix / 1000

# rego = ctem.util.read_linked_list_binary(path + '/rego_000001.dat', gridsize)
# melt = ctem.util.read_linked_list_binary(path + '/melt_000001.dat', gridsize)
# #age = ctem.util.read_linked_list_binary(path + '/age_000001.dat', gridsize, kind='SP')
# stack = ctem.util.read_unformatted_binary(path + '/stack_000001.dat', gridsize, kind='I4B')

rego = ctem.util.read_linked_list_binary(path + '/surface_rego.dat', gridsize)
melt = ctem.util.read_linked_list_binary(path + '/surface_melt.dat', gridsize)
#age = ctem.util.read_linked_list_binary(path + '/surface_age.dat', gridsize, kind='SP')
stack = ctem.util.read_unformatted_binary(path + '/surface_stacknum.dat', gridsize, kind='I4B')

regodepth = np.zeros((gridsize,gridsize))
for i in np.arange(gridsize):
    for j in np.arange(gridsize):
        regodepth[i,j] = np.sum(rego[i,j][2:]) / 1000

meltdepth = np.zeros((gridsize,gridsize))
for i in np.arange(gridsize):
    for j in np.arange(gridsize):
        meltdepth[i,j] = np.sum(melt[i,j][2:]) / stack[i,j]

#agearr = np.reshape(age[255,300], (60,196), order='F')
# agearr2 = np.reshape(age[140,150], (60,86), order='F')
# agearr3 = np.reshape(age[200,90], (60,7), order='F')

dpi = 72
x = np.arange(0,1000,200)

fig = plt.figure(figsize=(gridsize/dpi, gridsize/dpi), dpi=dpi)
ax = plt.axes([0,0,1,1])
rd = ax.imshow(regodepth, interpolation='nearest')#, vmin=0.0, vmax=10.0)
cbar = fig.colorbar(rd, label = 'Regolith Depth (km)')#, shrink=0.8)
ax.invert_yaxis()
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
xt = ax.get_xticks().tolist()
xt[1] = pixkm*x[1]
xt[2] = pixkm*x[2]
xt[3] = pixkm*x[3]
xt[4] = pixkm*x[4]
yt = ax.get_yticks().tolist()
yt[1] = pixkm*x[1]
yt[2] = pixkm*x[2]
yt[3] = pixkm*x[3]
yt[4] = pixkm*x[4]
ax.set_xticklabels(xt)
ax.set_yticklabels(yt)

fig2 = plt.figure(figsize=(gridsize/dpi, gridsize/dpi), dpi=dpi)
ax2 = plt.axes([0,0,1,1])
md = ax2.imshow(meltdepth, interpolation='nearest')#, vmin=0.0, vmax=0.5)
cbar = fig2.colorbar(md, label = 'Melt Fraction')#, shrink=0.8)
ax2.invert_yaxis()

#xax = np.arange(60)
#f = plt.figure()
#a = f.add_subplot(111)
#a.bar(xax, (agearr[:,0]/stack[255,300]))
# f, a = plt.subplots(2,2)
# a[0,0].bar(xax, agearr[:,73])#[0:60])
# a[0,0].set_xlabel('Age bin')
# a[0,0].set_ylabel('Amount of melt')
# a[0,1].bar(xax, agearr2[:,8])#[0:60])
# a[1,0].bar(xax, agearr2[:,9])#[0:60])
# a[1,1].bar(xax, agearr3[:,0])#[0:60])
plt.show()