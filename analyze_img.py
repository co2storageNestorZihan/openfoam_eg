#%%
import gzip
import numpy as np
import matplotlib.pyplot as plt
import porespy as ps


# Path to your .raw.gz file
file_path = "/Users/zihanren/active/openfoam_eg/large3d/Bentheimer_1000c_3p0035um.raw.gz"

# You need to know these parameters for your specific raw file
# For a CT scan like Bentheimer sample, these might be:
width = 1000  # Width in pixels
height = 1000  # Height in pixels
depth = 1000  # Number of slices/layers
dtype = np.uint8  # Data type (common types: uint8, uint16, float32)

# Read the compressed file
with gzip.open(file_path, 'rb') as f:
    # Read the decompressed content
    data = f.read()
    
    # Convert to numpy array with proper shape
    volume = np.frombuffer(data, dtype=dtype)
    print(volume.shape) 
    print(volume.dtype)
    print(volume.min())
    print(volume.max())
    
    # Reshape to 3D volume
    # Note: You may need to adjust order depending on how data is stored
    volume = volume.reshape((depth, height, width))
    volume = volume < 1
    # calculate the porosity
    porosity = ps.metrics.porosity(volume)
    print(f"Porosity: {porosity}")

    # slice subvolume to 200 cubic
    subvolume = volume[500:700, 500:700, 500:700]
    porosity = ps.metrics.porosity(subvolume)
    print(f"Porosity of subvolume: {porosity}")
    
    plt.imshow(subvolume[2], cmap='gray')
    plt.show()

    # save subvolume as numpy array
    # make True as 255 and False as 0
    subvolume = subvolume.astype(np.uint8) * 255
    # save this subvolume as npy file
    np.save('large3d/subvolume.npy', subvolume)

#%%
