from horizon_utils import HorizonMesh

segy_file = '../../../DDW data/DDW_VR.sgy'
horizon_point_cloud_files = [
    '../../../DDW data/DDW_H00_VR',
    '../../../DDW data/DDW_H20_VR',
    '../../../DDW data/DDW_H40_VR',
]

std_amplitude = HorizonMesh.get_std_amplitude_from_segy(segy_file)
print("Standard Deviation of Amplitude:", std_amplitude)

for file in horizon_point_cloud_files:
    horizon_mesh = HorizonMesh.HorizonMesh.from_text_file(file, crs='utm', sep=r'\s+')
    horizon_mesh.initialise_transformations(segy_file)
    horizon_mesh.convert_to_ilcl()
    horizon_mesh.rescale_to_segy(segy_file)
    horizon_mesh.swap_yz()
    output_file = file + '.obj'
    horizon_mesh.to_file(output_file, format='trimesh')
    print(horizon_mesh.vertices.min(axis=0), horizon_mesh.vertices.max(axis=0))
    print(f"Saved transformed mesh to {output_file}")
    
    