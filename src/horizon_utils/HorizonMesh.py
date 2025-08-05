import numpy as np
import trimesh
import segyio
from scipy.spatial import Delaunay
from scipy.interpolate import griddata
import pandas as pd

def derive_transformations(segy_file):
    """Derive transformation parameters from a SEGY file."""
    with segyio.open(segy_file, "r", ignore_geometry=True) as f:
        inlines = f.attributes(segyio.TraceField.INLINE_3D)[:]
        crosslines = f.attributes(segyio.TraceField.CROSSLINE_3D)[:]
        xs = f.attributes(segyio.TraceField.CDP_X)[:]/100
        ys = f.attributes(segyio.TraceField.CDP_Y)[:]/100

    x_utm = np.column_stack([xs, ys])
    x_ilcl = np.column_stack([inlines, crosslines])

    i0 = 0
    i1 = np.argmax(inlines)
    i2 = np.argmax(crosslines)

    assert inlines[i0] == inlines[i2]
    assert crosslines[i0] == crosslines[i1]

    stretch_matrix = np.array([
        [np.linalg.norm(x_utm[i1] - x_utm[i0]) / np.linalg.norm(x_ilcl[i1] - x_ilcl[i0]), 0],
        [0, np.linalg.norm(x_utm[i2] - x_utm[i0]) / np.linalg.norm(x_ilcl[i2] - x_ilcl[i0])]
    ])

    rotation_angle = np.arctan2((x_utm[i1,1] - x_utm[i0,1]), (x_utm[i1,0] - x_utm[i0,0]))

    rotation_matrix = np.array([
        [np.cos(rotation_angle), -np.sin(rotation_angle)],
        [np.sin(rotation_angle), np.cos(rotation_angle)]
    ])

    origin_ilcl = x_ilcl[i0]
    origin_utm = x_utm[i0]

    def transform_ix_to_utm(x_ilcl):
        """Transform (inline, crossline) to (easting, northing)."""
        return origin_utm[None,:] + (rotation_matrix @ stretch_matrix @ (x_ilcl - origin_ilcl).T).T
    
    def transform_utm_to_ix(x_utm):
        """Transform (easting, northing) to (inline, crossline)."""
        return origin_ilcl[None,:] + (np.linalg.inv(rotation_matrix @ stretch_matrix) @ (x_utm - origin_utm).T).T
    
    def transform_ix_to_dimensioned_ix(x_ilcl):
        """Transform (inline, crossline) to dimensioned inline/crossline."""
        return x_ilcl @ stretch_matrix.T

    return transform_ix_to_utm, transform_utm_to_ix, transform_ix_to_dimensioned_ix


class HorizonMesh:
    def __init__(self, vertices, faces, crs):
        self.vertices = vertices
        self.faces = faces
        self.crs = crs

        # Initialise transformation functions
        self.transform_ix_to_utm = None
        self.transform_utm_to_ix = None
        self.transform_ix_to_dimensioned_ix = None

    def initialise_transformations(self, segy_file):
        """Initialise the transformation functions for inline/crossline to UTM and vice versa."""
        self.transform_ix_to_utm, self.transform_utm_to_ix, self.transform_ix_to_dimensioned_ix = derive_transformations(segy_file)

    def convert_to_utm(self):
        """Transform mesh vertices from inline/crossline to UTM coordinates."""
        if self.crs == 'utm':
            print("Mesh is already in UTM coordinates.")
            return
        elif self.transform_ix_to_utm is None:
            raise ValueError("Transformation functions not initialized. Call initialise_transformations with a SEGY file first.")
        elif self.crs == 'inline_crossline':
            self.vertices[:,:2] = self.transform_ix_to_utm(self.vertices[:,:2])
            self.crs = 'utm'
        else:
            raise NotImplementedError("Conversion from current CRS to UTM is not implemented.")

    def convert_to_ilcl(self):
        """Transform mesh vertices from UTM coordinates to inline/crossline."""
        if self.crs == 'inline_crossline':
            print("Mesh is already in inline/crossline coordinates.")
            return
        elif self.transform_utm_to_ix is None:
            raise ValueError("Transformation functions not initialized. Call initialise_transformations with a SEGY file first.")
        elif self.crs == 'utm':
            self.vertices[:,:2] = self.transform_utm_to_ix(self.vertices[:,:2])
            self.crs = 'inline_crossline'
        else:
            raise NotImplementedError("Conversion from current CRS to inline/crossline is not implemented.")
    
    def convert_to_dimensioned_ilcl(self):
        """Transform mesh vertices to dimensioned inline/crossline coordinates."""
        if self.crs == 'dimensioned_inline_crossline':
            print("Mesh is already in dimensioned inline/crossline coordinates.")
            return
        elif self.transform_ix_to_dimensioned_ix is None:
            raise ValueError("Transformation functions not initialized. Call initialise_transformations with a SEGY file first.")
        elif self.crs == 'inline_crossline':
            self.vertices[:,:2] = self.transform_ix_to_dimensioned_ix(self.vertices[:,:2])
            self.crs = 'dimensioned_inline_crossline'
        elif self.crs == 'utm':
            self.convert_to_ilcl()
            self.vertices[:,:2] = self.transform_ix_to_dimensioned_ix(self.vertices[:,:2])
            self.crs = 'dimensioned_inline_crossline'
        else:
            raise NotImplementedError
    
    @classmethod
    def from_ts(cls, ts_path):
        """Load mesh from a TS file (e.g. generated by Petrel)

        Args:
            ts_path: Path to the TS file.

        Returns:
            HorizonMesh: An instance of HorizonMesh containing vertices and faces.
        """
        vertices = []
        faces = []

        with open(ts_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                if parts[0] in ('VRTX', 'PVRTX'):
                    vertices.append([float(parts[2]), float(parts[3]), float(parts[4])])
                elif parts[0] == 'TRGL':
                    faces.append([int(parts[1]), int(parts[2]), int(parts[3])])

        return cls(np.array(vertices), np.array(faces), crs='utm')

    @classmethod
    def from_text_file(cls, file_path, crs, **kwargs):
        """Load mesh from a text file.

        Args:
            file_path: Path to the text file containing vertices.
            **kwargs: Additional arguments for loading vertices.

        Returns:
            HorizonMesh: An instance of HorizonMesh.
        """
        vertices = pd.read_csv(file_path, **kwargs).values
        return cls.from_vertices(vertices, crs)
    
    @classmethod
    def from_DUG_ilclt(cls, filename):
        """Load mesh from a DUG ilclt file.

        Args:
            filename: Path to the DUG ilclt file.

        Returns:
            HorizonMesh: An instance of HorizonMesh.
        """
        return cls.from_text_file(filename, crs='inline_crossline', sep=r'\s+', comment='#')
    
    @classmethod
    def from_dugmsh(cls, filename):
        """Load mesh from a .dugmsh file.

        Args:
            filename: Path to the dugmsh file.

        Returns:
            HorizonMesh: An instance of HorizonMesh.
        """
        return cls.from_text_file(filename, crs='inline_crossline', sep=r'\s+', skiprows=14)
    
    @classmethod
    def from_DUG_dat(cls, filename):
        """Load mesh from a DUG dat file.

        Args:
            filename: Path to the DUG dat file.

        Returns:
            HorizonMesh: An instance of HorizonMesh.
        """
        return cls.from_text_file(filename, crs='inline_crossline', sep=r'\s+', usecols=[0, 1, 4])
    
    @classmethod
    def from_vertices(cls, vertices, crs):
        """Create a HorizonMesh instance from vertices.

        Args:
            vertices: A numpy array of shape (n, 3) representing the mesh vertices.
            crs: The coordinate reference system of the vertices.

        Returns:
            HorizonMesh: An instance of HorizonMesh.
        """
        tri = Delaunay(vertices[:,:2])
        faces = tri.simplices
        return cls(vertices, faces, crs=crs)
        
    def to_file(self, filename, format, **kwargs):
        """Save the HorizonMesh to a file.

        Args:
            filename: Path to the output file.
            format: The format to save the file in (e.g. 'trimesh', 'text').

        Raises:
            NotImplementedError: If the format is not supported.
        """
        if format == 'trimesh':
            m = trimesh.Trimesh(vertices=self.vertices, faces=self.faces, process=False)
            m.export(filename)
        elif format == 'text':
            with open(filename, 'w') as f:
                sep = kwargs.get('sep', '\t')
                vertices_df = pd.DataFrame(self.vertices, columns=['x', 'y', 'z'])
                vertices_df.to_csv(f, sep=sep, index=False, header=False)
        else:
            raise NotImplementedError

    def regrid(self, resolution):
        """Regrid the mesh to a new resolution.

        Args:
            resolution: The new resolution for the mesh.
        """
        new_x = np.arange(int(self.vertices[:, 0].min())+1, self.vertices[:, 0].max(), resolution)
        new_y = np.arange(int(self.vertices[:, 1].min())+1, self.vertices[:, 1].max(), resolution)
        new_x, new_y = np.meshgrid(new_x, new_y)

        # Regrid
        new_z = griddata((self.vertices[:,0], self.vertices[:,1]), self.vertices[:, 2],
                          (new_x, new_y), method='linear')
        
        non_nan_mask = ~np.isnan(new_z)
        new_z = new_z[non_nan_mask]
        new_x = new_x[non_nan_mask]
        new_y = new_y[non_nan_mask]

        new_vertices = np.column_stack([new_x.flatten(), new_y.flatten(), new_z.flatten()])

        # Re-mesh with new vertices
        new_faces = Delaunay(new_vertices[:, :2]).simplices
        
        self.vertices = new_vertices
        self.faces = new_faces

        return self
    
    def rescale(self, target_size=[1,1,1]):
        """Rescale the mesh to a new size.

        Args:
            target_size: A list of three values representing the target size in each dimension.
        """
        scale_factors = np.array(target_size) / np.array(self.vertices.max(axis=0) - self.vertices.min(axis=0))
        self.vertices *= scale_factors

        return self
    
    def reset_origin(self, new_origin=[0,0,0]):
        """Reset the origin of the mesh to a new point.

        Args:
            new_origin: A list of three values representing the new origin coordinates.
        """
        self.vertices -= self.vertices.mean(axis=0) - np.array(new_origin)

        return self
    
    def rescale_to_segy(self, segy_file):
        """Rescale all axes to 0-1 range based on seismic volume provided.

        Args:
            segy_file: Path to the SEGY file to derive the scaling factors.
        """
        if not self.crs == 'inline_crossline':
            self.convert_to_ilcl()
            print('Converted to inline/crossline coordinates prior to rescaling.')
        with segyio.open(segy_file, "r", ignore_geometry=True) as f:
            inlines = f.attributes(segyio.TraceField.INLINE_3D)[:]
            crosslines = f.attributes(segyio.TraceField.CROSSLINE_3D)[:]
            x_min = inlines.min()
            x_max = inlines.max()
            y_min = crosslines.min()
            y_max = crosslines.max()
                
            sample_interval = f.bin[segyio.BinField.Interval] / 1e6
            num_samples = f.bin[segyio.BinField.Samples]
            zs = np.arange(num_samples) * sample_interval
            z_min = zs.min()
            z_max = zs.max()
        
        self.vertices[:, 0] = (self.vertices[:, 0] - x_min) / (x_max - x_min)
        self.vertices[:, 1] = (self.vertices[:, 1] - y_min) / (y_max - y_min)
        self.vertices[:, 2] = (self.vertices[:, 2] - z_min) / (z_max - z_min)

        self.crs = 'scaled_to_segy'

        return self