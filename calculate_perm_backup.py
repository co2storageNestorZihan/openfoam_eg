#!/usr/bin/env python
# -*- coding: utf-8 -*-

import vtk
import numpy as np
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# Create output directory for figures if it doesn't exist
if not os.path.exists('figures'):
    os.makedirs('figures')

# Function to process a VTK file and calculate permeability
def process_vtk(vtkFile, time_step):
    print(f"\nProcessing time step {time_step}...")
    
    # load VTK data
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFile)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    
    # Calculate cell volumes
    cellVolumes = vtk.vtkCellSizeFilter()
    cellVolumes.SetInputData(data)
    cellVolumes.SetComputeVolume(True)
    cellVolumes.Update()
    volumeData = cellVolumes.GetOutput()
    Vcells = VN.vtk_to_numpy(volumeData.GetCellData().GetArray("Volume"))
    
    # extract velocity arrays
    U = VN.vtk_to_numpy(data.GetCellData().GetArray('U'))
    Umag = np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
    Ux = U[:,0]
    Uy = U[:,1]
    Uz = U[:,2]
    
    # extract cell centers for visualization
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(data)
    cell_centers.Update()
    centers = cell_centers.GetOutput()
    points = VN.vtk_to_numpy(centers.GetPoints().GetData())
    x_coords = points[:, 0]
    y_coords = points[:, 1]
    
    # calculate volume fluxes
    qx = []
    qy = []
    for c in range(len(Vcells)):
        q1 = Vcells[c] * Ux[c]
        q2 = Vcells[c] * Uy[c]
        qx.append(q1)
        qy.append(q2)
    
    # simulation parameters
    DP = 1                       # pressure drop [Pa]
    nu = 1e-06                   # kinematic viscosity [m²/s²]
    rho = 1000                   # density [kg/m³]
    mu = nu * rho                # dynamic viscosity [kg/(m*s)]
    dx = 0.001196                # model length x [m]
    dy = 0.001494                # model width y [m]
    dz = 1e-6                    # model thickness z [m]
    A = dy * dz
    V = dx * A
    
    # calculate Darcy velocity
    U_Darcy_x = np.sum(qx)/V
    
    kxx = mu * dx * U_Darcy_x/DP
    
    # Convert permeability to millidarcy (md)
    md_conversion = 1.01325e15  # 1 m² = 1.01325 × 10¹⁵ md
    kxx_md = kxx * md_conversion
    
    print(f'Time step {time_step} - Bulk permeability: {kxx:.3e} m² ({kxx_md:.3e} md)')
    
    # Create velocity plots
    plt.figure(figsize=(12, 10))
    
    # Plot 1: Velocity magnitude
    plt.subplot(221)
    scatter = plt.scatter(x_coords, y_coords, c=Umag, s=1, cmap='jet')
    plt.title(f'Time {time_step} - Velocity Magnitude')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(scatter, cax=cax)
    
    # Plot 2: X-component of velocity
    plt.subplot(222)
    scatter = plt.scatter(x_coords, y_coords, c=Ux, s=1, cmap='jet')
    plt.title(f'Time {time_step} - X-component of Velocity')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(scatter, cax=cax)
    
    # Plot 3: Y-component of velocity
    plt.subplot(223)
    scatter = plt.scatter(x_coords, y_coords, c=Uy, s=1, cmap='jet')
    plt.title(f'Time {time_step} - Y-component of Velocity')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(scatter, cax=cax)
    
    # Plot 4: Velocity vector field (downsampled for clarity)
    plt.subplot(224)
    # Downsample for clearer vector plot
    downsample = max(1, len(x_coords) // 2000)  # Adjust based on point density
    plt.quiver(x_coords[::downsample], y_coords[::downsample], 
               Ux[::downsample], Uy[::downsample], 
               Umag[::downsample], cmap='jet', scale=10)
    plt.title(f'Time {time_step} - Velocity Vector Field')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.colorbar()
    
    plt.tight_layout()
    plt.savefig(f'figures/velocity_maps_time_{time_step}.png', dpi=300)
    plt.close()
    
    return {
        'kxx': kxx,
        'kxx_md': kxx_md,
        'x_coords': x_coords,
        'y_coords': y_coords,
        'Umag': Umag,
        'Ux': Ux,
        'Uy': Uy
    }

# Process both time steps
results_0 = process_vtk('VTK/case_0.vtk', 0)
results_537 = process_vtk('VTK/case_537.vtk', 537)

# Create comparison plots
plt.figure(figsize=(15, 10))

# Plot Velocity Magnitude Comparison
plt.subplot(221)
plt.scatter(results_0['x_coords'], results_0['y_coords'], c=results_0['Umag'], s=1, cmap='jet', label='Time 0')
plt.title('Time 0 - Velocity Magnitude')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar()

plt.subplot(222)
plt.scatter(results_537['x_coords'], results_537['y_coords'], c=results_537['Umag'], s=1, cmap='jet', label='Time 537')
plt.title('Time 537 - Velocity Magnitude')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar()

# Calculate velocity difference for comparison
# Note: For simplicity, we assume the mesh coordinates are the same between time steps
# In a real application, you might need to interpolate between different meshes
plt.subplot(223)
# Create a histogram of velocity magnitudes
plt.hist(results_0['Umag'], bins=50, alpha=0.5, label='Time 0')
plt.hist(results_537['Umag'], bins=50, alpha=0.5, label='Time 537')
plt.title('Velocity Magnitude Distribution Comparison')
plt.xlabel('Velocity Magnitude [m/s]')
plt.ylabel('Frequency')
plt.legend()

plt.subplot(224)
# Show permeability comparison as a bar chart
plt.bar(['Time 0', 'Time 537'], [results_0['kxx'], results_537['kxx']])
plt.title('Permeability Comparison')
plt.ylabel('Permeability [m²]')
for i, v in enumerate([results_0['kxx'], results_537['kxx']]):
    plt.text(i, v, f'{v:.2e}', ha='center', va='bottom')

plt.tight_layout()
plt.savefig('figures/time_step_comparison.png', dpi=300)
plt.close()

print("\nPermeability Comparison:")
print(f"Time 0:   {results_0['kxx']:.3e} m² ({results_0['kxx_md']:.3e} md)")
print(f"Time 537: {results_537['kxx']:.3e} m² ({results_537['kxx_md']:.3e} md)")
print(f"Ratio (537/0): {results_537['kxx']/results_0['kxx']:.3f}")
print("\nFigures saved in the 'figures' directory.") 