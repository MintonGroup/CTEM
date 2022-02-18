import numpy as np
import os
import open3d
import ctem_io_readers
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
        ygrid, xgrid = np.mgrid[-s / 2:s / 2, -s / 2:s / 2] * pix
    
        xvals = xgrid.flatten()
        yvals = ygrid.flatten()
        zvals = self.dem.flatten()
        verts = np.array((xvals, yvals, zvals)).T
        faces = np.empty([2 * (s - 1) ** 2, 3], dtype=np.int64)
        for j in range(s - 1):
            for i in range(s - 1):
                idx = (s - 1) * j + i
                # faces[idx,:] = [ j*s+i, j*s+i+1, (j+1)*s+i+1, (j+1)*s+i ]
                faces[idx, :] = [j * s + i, j * s + i + 1, (j + 1) * s + i + 1]
                idx += (s - 1) ** 2
                faces[idx, :] = [(j + 1) * s + i + 1, (j + 1) * s + i, j * s + i ]
        self.mesh.vertices = open3d.utility.Vector3dVector(verts)
        self.mesh.triangles = open3d.utility.Vector3iVector(faces)
        self.mesh.compute_vertex_normals()
        
        return
    
    def visualize(self):
        vis = open3d.visualization.Visualizer()
        vis.create_window()
        vis.add_geometry(self.mesh)
        

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

    

    

