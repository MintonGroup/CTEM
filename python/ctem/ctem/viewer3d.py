import numpy as np
import open3d
import ctem

class Polysurface(ctem.Simulation):
    """A model of a self-gravitating small body"""
    def __init__(self):
        ctem.Simulation.__init__(self, isnew=False) # Initialize the data structures, but doesn't start a new run
        self.read_output()
        #used for Open3d module
        self.mesh = open3d.geometry.TriangleMesh()
        return
    
    def ctem2blockmesh(self):
        # Build mesh grid
        s = self.user['gridsize']
        pix = self.user['pix']
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
                    faces[fidx[2],:] = [i0, i3, i2-1]
                    faces[fidx[3],:] = [i0, i2-1, i1-1]
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
    
    def ctem2trimesh(self):
        # Build mesh grid
        s = self.user['gridsize']
        pix = self.user['pix']
        ygrid, xgrid = np.mgrid[-s/2:s/2, -s/2:s/2] * pix

        xvals = xgrid.flatten()
        yvals = ygrid.flatten()
        zvals = self.surface_dem.flatten()
        verts = np.array((xvals, yvals, zvals)).T
        faces = np.full([2 * (s-1)**2, 3], -1, dtype=np.int64)
        for j in range(s - 1):
            for i in range(s - 1):
                idx = (s - 1) * j + i
                faces[idx, :] = [j * s + i, j * s + i + 1, (j + 1) * s + i + 1]
                idx += (s - 1) ** 2
                faces[idx, :] = [(j + 1) * s + i + 1, (j + 1) * s + i, j * s + i ]
        self.mesh.vertices = open3d.utility.Vector3dVector(verts)
        self.mesh.triangles = open3d.utility.Vector3iVector(faces)
        self.mesh.compute_vertex_normals()
        self.mesh.compute_triangle_normals()

        return
    
    def visualize(self):
        vis = open3d.visualization.Visualizer()
        vis.create_window()
        vis.add_geometry(self.mesh)
        opt = vis.get_render_option()
        opt.background_color = np.asarray([0, 0, 0])

        self.mesh.paint_uniform_color([0.5, 0.5, 0.5])
        vis.run()
        vis.destroy_window()


if __name__ == '__main__':
    sim = Polysurface()
    sim.ctem2trimesh()
    sim.visualize()

    

