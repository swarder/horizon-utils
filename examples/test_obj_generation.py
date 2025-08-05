from horizon_utils import HorizonMesh

segy_file = 'DDW_subcube_1.sgy'

mesh = HorizonMesh.HorizonMesh.from_DUG_ilclt('H53_subcube_mesh.ilclt',)
print(mesh.crs)
mesh.to_file('H53_subcube_mesh_ilclt.obj', format='trimesh')

mesh = HorizonMesh.HorizonMesh.from_dugmsh('H53_subcube_mesh.dugmsh')
print(mesh.crs)
mesh.to_file('H53_subcube_mesh_dugmsh.obj', format='trimesh')

mesh = HorizonMesh.HorizonMesh.from_DUG_dat('H53_subcube_mesh.dat')
print(mesh.crs)
mesh.to_file('H53_subcube_mesh_dat.obj', format='trimesh')

mesh = HorizonMesh.HorizonMesh.from_ts('H53_subcube_mesh.ts')
mesh.initialise_transformations(segy_file)
mesh.convert_to_dimensioned_ilcl()
print(mesh.crs)
mesh.to_file('H53_subcube_mesh_dimensioned_ilcl.obj', format='trimesh')

print(mesh.vertices.mean(axis=0))


mesh = HorizonMesh.HorizonMesh.from_ts('H53_subcube_mesh.ts')
mesh.initialise_transformations(segy_file)
mesh.rescale_to_segy(segy_file)
mesh.to_file('H53_subcube_mesh_rescaled_to_volume.obj', format='trimesh')