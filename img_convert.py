import numpy as np
import porespy as ps
from paraview.simple import *

# 1. Load 3D numpy array from .npy file
npy_path = 'large3d/subvolume.npy'
arr_3d = np.load(npy_path)

# Check the shape of the array
print(f"Array shape: {arr_3d.shape}")
print(f"Array min value: {np.min(arr_3d)}, max value: {np.max(arr_3d)}")

# 2. Save as vti (VTK's image format) using porespy
# No need to stack since the array is already 3D
ps.io.to_vtk(arr_3d, 'geometry/porous_model')


def write_stl(vti_file, stl_file):
    # 1. load vti file
    data            = OpenDataFile('geometry/%s.vti'%vti_file)
    # 2. clip at some intermediate value (we have 0 and 255 as solid and pores)
    # Invert=0 to keep the solid part (value 0) and remove pores (value 255)
    clip1           = Clip(data, ClipType = 'Scalar', Scalars = ['CELLS', 'im'], Value = 127.5, Invert = 0)
    # 3. make a surface of the remaining solid parts
    extractSurface1 = ExtractSurface(clip1)
    # 4. and triangulate it for stl export
    triangulate1    = Triangulate(extractSurface1)

    # 5. finally save it as an stl file
    SaveData(stl_file, proxy = triangulate1)

# main part
vti_file = 'porous_model'      # input .vti file
stl_file = 'porous_model.stl'  # output .stl file
# call function
write_stl(vti_file, stl_file)