import numpy as np
from ctem import ctem_viewer_3d

surf = ctem_viewer_3d.surf()
surf.io_read_ply("216kleopatra-mesh.ply")
v = np.asarray(surf.mesh)
t = np.asarray(surf.mesh)
surf.render()