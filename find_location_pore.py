import numpy as np
# You might need to install scipy: pip install scipy
from scipy.ndimage import distance_transform_edt

# 1. Load your original 3D numpy array
npy_path = 'large3d/subvolume.npy'
arr_3d = np.load(npy_path)

# 2. Create a boolean array representing the pore space (True where pores exist)
#    Assuming pores are 255
is_pore = (arr_3d == 255)

# 3. Calculate the distance transform
#    For each True voxel (pore), this finds the distance to the nearest False voxel (solid)
distance_map = distance_transform_edt(is_pore)

# 4. Find the indices of the pore voxel with the maximum distance
#    This is the voxel "deepest" inside the pore space
max_dist_indices = np.unravel_index(np.argmax(distance_map), distance_map.shape)
i, j, k = max_dist_indices

# --- VITAL ASSUMPTION ---
# Assuming porespy saved the VTI with origin (0,0,0) and spacing (1,1,1)
# PLEASE VERIFY THIS: Open geometry/porous_model.vti in ParaView,
# check the 'Information' panel for 'Origin' and 'Spacing'.
# If they are different, adjust the calculation below.
origin = (0.0, 0.0, 0.0)
spacing = (1.0, 1.0, 1.0)
# --- END ASSUMPTION ---

# 5. Calculate the physical coordinates of the center of this voxel
#    This is the point snappyHexMesh needs
x = origin[0] + (i + 0.5) * spacing[0]
y = origin[1] + (j + 0.5) * spacing[1]
z = origin[2] + (k + 0.5) * spacing[2]

print(f"Original array shape: {arr_3d.shape}")
print(f"Pore voxel indices with max distance from solid: ({i}, {j}, {k})")
print(f"Distance value: {distance_map[i, j, k]}")
print(f"Recommended locationInMesh: ({x} {y} {z})")

# For direct copy-pasting into snappyHexMeshDict:
print(f"\nlocationInMesh ({x} {y} {z});")