import numpy as np
from ctem import ctem_viewer_3d

surf = ctem_viewer_3d.surf()
surf.io_read_ply("216kleopatra-mesh.ply")
v = np.asarray(surf.o3d_mesh)
t = np.asarray(surf.o3d_mesh)
surf.render()