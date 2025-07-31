from horizon_utils import HorizonMesh

mesh = HorizonMesh.HorizonMesh.from_DUG_ilclt('H53_subcube_mesh.ilclt')
print(mesh.vertices)

mesh = HorizonMesh.HorizonMesh.from_dugmsh('H53_subcube_mesh.dugmsh')
print(mesh.vertices)

mesh = HorizonMesh.HorizonMesh.from_DUG_dat('H53_subcube_mesh.dat')
print(mesh.vertices)