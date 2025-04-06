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
    
    # extract pressure data
    p = VN.vtk_to_numpy(data.GetCellData().GetArray('p'))
    
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
    
    # Determine inlet boundary cells (assuming inlet is at minimum x)
    inlet_cells = x_coords < min(x_coords) + 0.0001
    inlet_pressure = np.mean(p[inlet_cells])
    inlet_velocity = np.mean(Ux[inlet_cells])
    
    # Print boundary conditions for time step 0
    if time_step == 0:
        print("\nInitial Boundary Conditions:")
        print(f"Inlet Pressure: {inlet_pressure:.6f} Pa")
        print(f"Inlet Velocity: {inlet_velocity:.6f} m/s")
    
    return {
        'kxx': kxx,
        'kxx_md': kxx_md,
        'x_coords': x_coords,
        'y_coords': y_coords,
        'Umag': Umag,
        'Ux': Ux,
        'Uy': Uy,
        'p': p
    }

# Process both time steps
results_0 = process_vtk('VTK/case_0.vtk', 0)
results_537 = process_vtk('VTK/case_537.vtk', 537)

# Figure 1: Velocity field spatial distribution for time 0 and time 537
plt.figure(figsize=(15, 7))

# Time 0 velocity plot
plt.subplot(121)
# Downsample for clearer vector plot
downsample = max(1, len(results_0['x_coords']) // 2000)
plt.quiver(results_0['x_coords'][::downsample], results_0['y_coords'][::downsample], 
           results_0['Ux'][::downsample], results_0['Uy'][::downsample], 
           results_0['Umag'][::downsample], cmap='jet', scale=10)
plt.title('Time 0 - Velocity Field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(label='Velocity Magnitude [m/s]')

# Time 537 velocity plot
plt.subplot(122)
downsample = max(1, len(results_537['x_coords']) // 2000)
plt.quiver(results_537['x_coords'][::downsample], results_537['y_coords'][::downsample], 
           results_537['Ux'][::downsample], results_537['Uy'][::downsample], 
           results_537['Umag'][::downsample], cmap='jet', scale=10)
plt.title('Time 537 - Velocity Field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(label='Velocity Magnitude [m/s]')

plt.tight_layout()
plt.savefig('figures/figure1_velocity_field.png', dpi=300)
plt.close()

# Figure 2: Pressure field spatial distribution for time 0 and time 537
plt.figure(figsize=(15, 7))

# Time 0 pressure plot
plt.subplot(121)
scatter = plt.scatter(results_0['x_coords'], results_0['y_coords'], c=results_0['p'], s=1, cmap='viridis')
plt.title('Time 0 - Pressure Field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(scatter, label='Pressure [Pa]')

# Time 537 pressure plot
plt.subplot(122)
scatter = plt.scatter(results_537['x_coords'], results_537['y_coords'], c=results_537['p'], s=1, cmap='viridis')
plt.title('Time 537 - Pressure Field')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(scatter, label='Pressure [Pa]')

plt.tight_layout()
plt.savefig('figures/figure2_pressure_field.png', dpi=300)
plt.close()

# Figure 3: Histogram plots of velocity at time 0 and time 537
plt.figure(figsize=(15, 7))

# Time 0 velocity histogram
plt.subplot(121)
plt.hist(results_0['Umag'], bins=50, color='blue')
plt.title('Time 0 - Velocity Magnitude Distribution')
plt.xlabel('Velocity Magnitude [m/s]')
plt.ylabel('Frequency')

# Time 537 velocity histogram
plt.subplot(122)
plt.hist(results_537['Umag'], bins=50, color='red')
plt.title('Time 537 - Velocity Magnitude Distribution')
plt.xlabel('Velocity Magnitude [m/s]')
plt.ylabel('Frequency')

plt.tight_layout()
plt.savefig('figures/figure3_velocity_histogram.png', dpi=300)
plt.close()

# Figure 4: Histogram plots of pressure at time 0 and time 537
plt.figure(figsize=(15, 7))

# Time 0 pressure histogram
plt.subplot(121)
plt.hist(results_0['p'], bins=50, color='blue')
plt.title('Time 0 - Pressure Distribution')
plt.xlabel('Pressure [Pa]')
plt.ylabel('Frequency')

# Time 537 pressure histogram
plt.subplot(122)
plt.hist(results_537['p'], bins=50, color='red')
plt.title('Time 537 - Pressure Distribution')
plt.xlabel('Pressure [Pa]')
plt.ylabel('Frequency')

plt.tight_layout()
plt.savefig('figures/figure4_pressure_histogram.png', dpi=300)
plt.close()

# Figure 5: Enhanced spatial velocity distribution showing the framework
plt.figure(figsize=(15, 12))

# Time 0 - Detailed spatial velocity distribution with framework
plt.subplot(221)
# Background scatter plot to show velocity magnitude (represents framework)
scatter = plt.scatter(results_0['x_coords'], results_0['y_coords'], 
                     c=results_0['Umag'], s=2, cmap='viridis', alpha=0.8)
plt.title('Time 0 - Velocity Magnitude Distribution')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(scatter, label='Velocity Magnitude [m/s]')

# Time 0 - Vector plot overlay
plt.subplot(222)
# Use tricontourf for better visualization of the flow field
from matplotlib.tri import Triangulation
triang = Triangulation(results_0['x_coords'], results_0['y_coords'])
contour = plt.tricontourf(triang, results_0['Umag'], 
                         levels=20, cmap='viridis', alpha=0.9)
# Add velocity vectors on top
downsample = max(1, len(results_0['x_coords']) // 1000)
plt.quiver(results_0['x_coords'][::downsample], results_0['y_coords'][::downsample], 
           results_0['Ux'][::downsample], results_0['Uy'][::downsample], 
           color='white', scale=15, width=0.002)
plt.title('Time 0 - Flow Field with Vectors')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(contour, label='Velocity Magnitude [m/s]')

# Time 537 - Detailed spatial velocity distribution with framework
plt.subplot(223)
# Background scatter plot to show velocity magnitude (represents framework)
scatter = plt.scatter(results_537['x_coords'], results_537['y_coords'], 
                     c=results_537['Umag'], s=2, cmap='viridis', alpha=0.8)
plt.title('Time 537 - Velocity Magnitude Distribution')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(scatter, label='Velocity Magnitude [m/s]')

# Time 537 - Vector plot overlay
plt.subplot(224)
# Use tricontourf for better visualization of the flow field
triang = Triangulation(results_537['x_coords'], results_537['y_coords'])
contour = plt.tricontourf(triang, results_537['Umag'], 
                         levels=20, cmap='viridis', alpha=0.9)
# Add velocity vectors on top
downsample = max(1, len(results_537['x_coords']) // 1000)
plt.quiver(results_537['x_coords'][::downsample], results_537['y_coords'][::downsample], 
           results_537['Ux'][::downsample], results_537['Uy'][::downsample], 
           color='white', scale=15, width=0.002)
plt.title('Time 537 - Flow Field with Vectors')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
plt.colorbar(contour, label='Velocity Magnitude [m/s]')

plt.tight_layout()
plt.savefig('figures/figure5_spatial_velocity_framework.png', dpi=300)
plt.close()

print("\nPermeability Comparison:")
print(f"Time 0:   {results_0['kxx']:.3e} m² ({results_0['kxx_md']:.3e} md)")
print(f"Time 537: {results_537['kxx']:.3e} m² ({results_537['kxx_md']:.3e} md)")
print(f"Ratio (537/0): {results_537['kxx']/results_0['kxx']:.3f}")
print("\nFigures saved in the 'figures' directory.")
