#!/usr/bin/env python
# -*- coding: utf-8 -*-

import vtk
import numpy as np
from vtk.util import numpy_support as VN
# import matplotlib.pyplot as plt # Removed
# from mpl_toolkits.axes_grid1 import make_axes_locatable # Removed
import os

# Create output directory for figures if it doesn't exist
# if not os.path.exists('figures'): # Removed plotting
#     os.makedirs('figures') # Removed plotting

# Function to process a VTK file and calculate permeability
def process_vtk(vtkFile, time_step):
    print(f"\nProcessing time step {time_step}...")

    # Check if file exists
    if not os.path.exists(vtkFile):
        print(f"Error: VTK file not found at {vtkFile}")
        return None

    # load VTK data
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFile)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()

    if data is None or data.GetNumberOfCells() == 0:
        print(f"Error: Failed to read data or no cells found in {vtkFile}")
        return None

    # Check if arrays exist
    if data.GetCellData().GetArray('U') is None:
        print(f"Error: 'U' array not found in cell data for {vtkFile}")
        return None
    if data.GetCellData().GetArray('p') is None:
        print(f"Error: 'p' array not found in cell data for {vtkFile}")
        return None

    # Calculate cell volumes
    cellVolumes = vtk.vtkCellSizeFilter()
    cellVolumes.SetInputData(data)
    cellVolumes.SetComputeVolume(True)
    cellVolumes.Update()
    volumeData = cellVolumes.GetOutput()
    Vcells = VN.vtk_to_numpy(volumeData.GetCellData().GetArray("Volume"))

    # extract velocity arrays
    U = VN.vtk_to_numpy(data.GetCellData().GetArray('U'))
    # Umag = np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2) # Not needed for calculation
    Ux = U[:,0]
    # Uy = U[:,1] # Not needed for calculation
    # Uz = U[:,2] # Not needed for calculation

    # extract pressure data
    p = VN.vtk_to_numpy(data.GetCellData().GetArray('p'))

    # extract cell centers for boundary check (optional, can be removed if not needed)
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(data)
    cell_centers.Update()
    centers = cell_centers.GetOutput()
    points = VN.vtk_to_numpy(centers.GetPoints().GetData())
    x_coords = points[:, 0]
    # y_coords = points[:, 1] # Not needed for calculation

    # calculate volume fluxes (only x-direction needed for kxx)
    qx_total = np.sum(Vcells * Ux)

    # simulation parameters
    DP = 1.0                     # pressure drop [Pa] (Confirmed)
    nu = 1e-06                   # kinematic viscosity [m²/s] (Confirmed)
    rho = 1000.0                 # density [kg/m³] (Assumed standard)
    mu = nu * rho                # dynamic viscosity [kg/(m*s)]

    # --- Corrected Physical Dimensions (based on voxel resolution) ---
    voxel_resolution = 3.0035e-6 # meters/voxel
    Nx = 200.0 # Number of voxels in x (adjust if different)
    Ny = 200.0 # Number of voxels in y (adjust if different)
    Nz = 200.0 # Number of voxels in z (adjust if different)

    dx = Nx * voxel_resolution   # model length x [m]
    dy = Ny * voxel_resolution   # model width y [m]
    dz = Nz * voxel_resolution   # model thickness z [m]
    print(f"Using Dimensions: dx={dx:.3e}m, dy={dy:.3e}m, dz={dz:.3e}m")

    # --- Cross-sectional area perpendicular to x-flow ---
    A = dy * dz
    # V = dx * A # Total Volume - not needed for Darcy velocity calculation

    # calculate Darcy velocity (Total Flux / Cross-Sectional Area)
    U_Darcy_x = qx_total / A

    # Calculate permeability using Darcy's Law: U = -(k/mu)*(dP/dx)
    # Rearranged: k = -U * mu * (dx/dP)
    # Note: We use absolute DP and U_Darcy_x, sign indicates flow direction
    kxx = mu * dx * U_Darcy_x / DP

    # Convert permeability to millidarcy (md)
    # 1 Darcy ≈ 0.9869233 × 10⁻¹² m² => 1 m² ≈ 1.00 / (0.9869233e-12) Darcy
    # 1 Darcy = 1000 mD
    # 1 m² ≈ (1.00 / 0.9869233e-12) * 1000 mD ≈ 1.01325e15 mD
    md_conversion = 1.0 / (0.9869233e-12 * 1e-3) # More precise calculation
    kxx_md = kxx * md_conversion

    print(f'Total X-Flux (Qx): {qx_total:.3e} m³/s')
    print(f'Cross-Sectional Area (A): {A:.3e} m²')
    print(f'Darcy Velocity (Ux): {U_Darcy_x:.3e} m/s')
    print(f'Time step {time_step} - Bulk permeability (kxx): {kxx:.3e} m² ({kxx_md:.3e} md)')

    # Determine inlet boundary cells (assuming inlet is at minimum x)
    # Be careful with floating point comparisons
    min_x = np.min(x_coords)
    # Use a small tolerance based on expected cell size (dx/200 for blockMesh)
    tolerance = dx / 200.0 * 0.5 # Tolerance = half a cell width
    inlet_cells = x_coords < (min_x + tolerance)

    if np.any(inlet_cells):
        inlet_pressure = np.mean(p[inlet_cells])
        inlet_velocity_x = np.mean(Ux[inlet_cells])
        # Print boundary conditions only for the first time step processed
        # if time_step == 0: # Let's print for both for comparison
        print(f"\nBoundary Check (Time Step {time_step}):")
        print(f"Approx Min X Coordinate: {min_x:.4f} m")
        print(f"Number of detected inlet cells: {np.sum(inlet_cells)}")
        print(f"Average Inlet Pressure: {inlet_pressure:.6f} Pa")
        print(f"Average Inlet X-Velocity: {inlet_velocity_x:.6f} m/s")
    else:
        print(f"Warning: Could not detect inlet cells near min x = {min_x:.4f} for time step {time_step}")


    return {
        'kxx': kxx,
        'kxx_md': kxx_md
        # Removed fields not needed for calculation
        # 'x_coords': x_coords,
        # 'y_coords': y_coords,
        # 'Umag': Umag,
        # 'Ux': Ux,
        # 'Uy': Uy,
        # 'p': p
    }

# Process both time steps
results_0 = process_vtk('VTK/case_0.vtk', 0)
results_537 = process_vtk('VTK/case_537.vtk', 537)

# --- Removed All Plotting Sections ---

# Final Permeability Comparison
print("\n" + "="*30)
print("Permeability Comparison:")
print("="*30)
if results_0 is not None:
    print(f"Time 0:   {results_0['kxx']:.3e} m² ({results_0['kxx_md']:.3e} md)")
else:
    print("Time 0:   Calculation failed.")

if results_537 is not None:
    print(f"Time 537: {results_537['kxx']:.3e} m² ({results_537['kxx_md']:.3e} md)")
else:
    print("Time 537: Calculation failed.")

if results_0 is not None and results_537 is not None and results_0['kxx'] != 0:
    print(f"Ratio (537/0): {results_537['kxx']/results_0['kxx']:.3f}")
elif results_0 is not None and results_0['kxx'] == 0:
     print("Ratio (537/0): Cannot calculate ratio (kxx at time 0 is zero).")

print("\nScript finished.")
