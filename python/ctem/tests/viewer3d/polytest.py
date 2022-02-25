import numpy as np
from ctem import ctem_viewer_3d
import open3d
surf = ctem_viewer_3d.polysurface()
#surf.io_read_ply("216kleopatra-mesh.ply")
surf.mesh = open3d.geometry.TriangleMesh.create_box()
v = np.asarray(surf.mesh.vertices)
t = np.asarray(surf.mesh.triangles)
surf.visualize()