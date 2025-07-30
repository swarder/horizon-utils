from horizon_utils import HorizonMesh
import numpy as np

# Define input/output paths
ts_file = 'H53_subcube_mesh.ts'
segy_file = 'DDW_subcube_1.sgy'
mesh_output_file = 'H53_subcube_mesh.obj'

# Load mesh from TS file
mesh = HorizonMesh.HorizonMesh.from_ts(ts_file)
print(mesh.crs)
print("Original vertices:\n", mesh.vertices)
orig_vertices = mesh.vertices.copy()

# Convert from UTM to inline/crossline using information from SEGY file
mesh.convert_from_utm(segy_file)
print(mesh.crs)
print("Vertices after converting from UTM:\n", mesh.vertices)

# Save the mesh to an OBJ file, in inline/crossline coordinates
mesh.save(mesh_output_file)

# Convert back to UTM coordinates to test transformation
mesh.convert_to_utm(segy_file)
print(mesh.crs)
print("Vertices after converting back to UTM:\n", mesh.vertices)

print('Recovered original vertices?')
print(np.allclose(orig_vertices, mesh.vertices))