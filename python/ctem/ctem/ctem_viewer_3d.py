import numpy as np
import os
import open3d
from ctem import ctem_io_readers
import os


class polysurface:
    """A model of a self-gravitating small body"""
    def __init__(self):
        self.polyfile = ""
        #used for Open3d module
        self.mesh = open3d.geometry.TriangleMesh()
        return
    
    def ctem2mesh(self, surface_file):
    
        # Create and initialize data dictionaries for parameters and options from CTEM.in
        notset = '-NOTSET-'
        currentdir = os.getcwd() + os.sep
    
        parameters = {'restart': notset,
                      'runtype': notset,
                      'popupconsole': notset,
                      'saveshaded': notset,
                      'saverego': notset,
                      'savepres': notset,
                      'savetruelist': notset,
                      'seedn': 1,
                      'totalimpacts': 0,
                      'ncount': 0,
                      'curyear': 0.0,
                      'fracdone': 1.0,
                      'masstot': 0.0,
                      'interval': 0.0,
                      'numintervals': 0,
                      'pix': -1.0,
                      'gridsize': -1,
                      'seed': 0,
                      'maxcrat': 1.0,
                      'shadedminhdefault': 1,
                      'shadedmaxhdefault': 1,
                      'shadedminh': 0.0,
                      'shadedmaxh': 0.0,
                      'workingdir': currentdir,
                      'ctemfile': 'ctem.in',
                      'datfile': 'ctem.dat',
                      'impfile': notset,
                      'sfdcompare': notset,
                      'sfdfile': notset}

        # Read ctem.in to initialize parameter values based on user input
        ctem_io_readers.read_ctemin(parameters, notset)
        # Read surface dem(shaded relief) and ejecta data files
        dem_file = parameters['workingdir'] + surface_file
        self.dem = ctem_io_readers.read_unformatted_binary(dem_file, parameters['gridsize'])
        # Build mesh grid
        s = parameters['gridsize']
        pix = parameters['pix']
        ygrid, xgrid = np.mgrid[-s/2:s/2, -s/2:s/2] * pix
        y0, x0 = ygrid - pix/2, xgrid - pix/2
        y1, x1 = ygrid - pix/2, xgrid + pix/2
        y2, x2 = ygrid + pix/2, xgrid + pix/2
        y3, x3 = ygrid + pix/2, xgrid - pix/2

    
        xvals = np.append(
                np.append(
                np.append(x0.flatten(),
                          x1.flatten()),
                          x2.flatten()),
                          x3.flatten())
        yvals = np.append(
                np.append(
                np.append(y0.flatten(),
                          y1.flatten()),
                          y2.flatten()),
                          y3.flatten())
        zvals = np.append(
                np.append(
                np.append(self.dem.flatten(),
                          self.dem.flatten()),
                          self.dem.flatten()),
                          self.dem.flatten())
        verts = np.array((xvals, yvals, zvals)).T
        nface_triangles = 10
        faces = np.full([nface_triangles*s**2, 3], -1, dtype=np.int64)
        for j in range(s):
            for i in range(s):
                i0 = s*j + i
                i1 = i0 + s**2
                i2 = i1 + s**2
                i3 = i2 + s**2

                fidx = np.arange(nface_triangles*i0,nface_triangles*(i0+1))
                # Make the two top triangles
                faces[fidx[0],:] = [i0, i1, i2]
                faces[fidx[1],:] = [i0, i2, i3]
                # Make the two west side triangles
                if i > 0:
                    faces[fidx[0],:] = [i0, i3, i3-1]
                    faces[fidx[1],:] = [i0, i3-1, i0-1]
                # Make the two south side triangles
                if j > 0:
                    faces[fidx[4],:] = [i1, i0, i3-s ]
                    faces[fidx[5],:] = [i1, i3-s, i2-s]
                # Make the two east side triangles
                if i < (s - 1):
                    faces[fidx[6],:] = [i2, i1, i0+1]
                    faces[fidx[7],:] = [i2, i0+1, i3+1]
                #Make the two north side triangles
                if j < (s -1):
                    faces[fidx[8],:] = [i3, i2, i1+s]
                    faces[fidx[9],:] = [i3, i1+s, i0+s]
        nz = faces[:,0] != -1
        f2 = faces[nz,:]
        self.mesh.vertices = open3d.utility.Vector3dVector(verts)
        self.mesh.triangles = open3d.utility.Vector3iVector(f2)
        self.mesh.compute_vertex_normals()
        self.mesh.compute_triangle_normals()

        return
    
    def visualize(self):
        vis = open3d.visualization.Visualizer()
        vis.create_window()
        vis.add_geometry(self.mesh)
        opt = vis.get_render_option()
        opt.background_color = np.asarray([0, 0, 0])

        # zmax = np.amax(surf.dem)
        # zmin = np.amin(surf.dem)
        # cnorm = (surf.dem.flatten() - zmin) / (zmax - zmin)
        # cval = plt.cm.terrain(cnorm)[:, :3]

        self.mesh.paint_uniform_color([0.5, 0.5, 0.5])
        vis.run()
        vis.destroy_window()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    surf = polysurface()
    surf.ctem2mesh(surface_file="surface_dem.dat")
    surf.visualize()
    #open3d.visualization.draw_geometries_with_editing([surf.mesh])

    

    

