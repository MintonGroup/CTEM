import numpy as np
import open3d as o3d
import trimesh

class surf:
    """A model of a self-gravitating small body"""
    def __init__(self):
        self.polyfile = ""
        self.mesh = trimesh.Trimesh()
        #used for Open3d module
        self.o3d_mesh = o3d.geometry.TriangleMesh()
        self.nface = 0
        self.nvert = 0
        return
    
    def io_read_ply(self,filename):
        """Reads in PLY file into Open3d triangular mesh object and computes normals"""
        self.polyfile = filename
        self.mesh = trimesh.load(f"{self.polyfile}")
        trimesh.repair.fix_normals(self.mesh)
        self.mesh.vertex_colors = np.full((self.nvert,4),1.0)
        self.nvert = len(self.mesh.vertices)
        self.nface = len(self.mesh.faces)
        self.convert_trimesh_to_o3d()
        return
    
    def render(self):
        """Renders output"""
        self.convert_trimesh_to_o3d()
        o3d.visualization.draw_geometries([self.o3d_mesh])
        return
    
    def convert_trimesh_to_o3d(self):
        """Update Open3D mesh with new vertices and faces, and recompute all normals"""
        self.o3d_mesh.vertices  = o3d.utility.Vector3dVector(self.mesh.vertices)
        self.o3d_mesh.triangles = o3d.utility.Vector3iVector(self.mesh.faces)
        self.o3d_mesh.compute_vertex_normals()
        self.o3d_mesh.compute_triangle_normals()
        self.o3d_mesh.vertex_colors = o3d.utility.Vector3dVector(self.mesh.vertex_colors[:,0:3])
        return

if __name__ == '__main__':
    import os
    import ctem_io_readers
    import open3d

    #Create and initialize data dictionaries for parameters and options from CTEM.in
    notset = '-NOTSET-'
    currentdir = os.getcwd() + os.sep

    parameters={'restart': notset,
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

    #Read ctem.in to initialize parameter values based on user input
    ctem_io_readers.read_ctemin(parameters,notset)
    #Read surface dem(shaded relief) and ejecta data files
    dem_file = parameters['workingdir'] + 'surface_dem.dat'
    surface_dem = ctem_io_readers.read_unformatted_binary(dem_file, parameters['gridsize'])
    # Build mesh grid
    gridsize = parameters['gridsize']
    pix = parameters['pix']
    ygrid, xgrid = np.mgrid[-gridsize / 2:gridsize / 2,
                             -gridsize / 2:gridsize / 2] * pix

    xvals = xgrid.flatten()
    yvals = ygrid.flatten()
    zvals = surface_dem.flatten()
    dem_vec = np.array((xvals, yvals, zvals*10)).T
    pcd = open3d.geometry.PointCloud()
    pcd.points = open3d.utility.Vector3dVector(dem_vec)

    o3d.visualization.draw_geometries([pcd])
    

