import numpy as np
from PIL import Image
import porespy as ps

impath = 'geometry/'         # image path
imname = 'porousModel.png'   # file name

# 1. use PIL to load the image
# convert("L") converts RGB into greyscale, L = R * 299/1000 + G * 587/1000 + B * 114/1000
image = Image.open('%s/%s'%(impath,imname)).convert("L")

# 2. convert image into numpy array with "normal" integer values
arr = np.asarray(image, dtype=int)

# check the shape of the array
print(arr.shape)

# 3. save as vti (VTK's image format)
# 3.1 make 3D by repeating in new axis=2 dimension (openfoam and porespy want 3D)
arr_stacked = np.stack((arr,arr), axis=2)

# 3.2 save as vti using porespy
ps.io.to_vtk(arr_stacked, 'geometry/porous_model')



from paraview.simple import *

def write_stl(vti_file, stl_file):
    # 1. load vti file
    data            = OpenDataFile('geometry/%s.vti'%vti_file)
    # 2. clip at some intermediate value (we have 0 and 255 as pores and grains)
    clip1           = Clip(data, ClipType = 'Scalar', Scalars = ['CELLS', 'im'], Value = 127.5, Invert = 1)
    # 3. make a surface of the remaining grains
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