#!/usr/bin/env python
# -*- coding: utf-8 -*-

#%%
import vtk
import numpy as np
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# Create output directory for 3D figures if it doesn't exist
if not os.path.exists('figures_3d'):
    os.makedirs('figures_3d')

def process_vtk_3d(vtkFile, time_step, include_solids=True):
    print(f"\nProcessing 3D data for time step {time_step}...")
    
    # load fluid VTK data
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
    Umag = np.sqrt(U[:,0]**2 + U[:,1]**2 + U[:,2]**2)
    
    # extract pressure data
    p = VN.vtk_to_numpy(data.GetCellData().GetArray('p'))
    
    # extract cell centers
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(data)
    cell_centers.Update()
    centers = cell_centers.GetOutput()
    fluid_points = VN.vtk_to_numpy(centers.GetPoints().GetData())
    
    # If requested, load the solid structure data
    solid_points = None
    if include_solids:
        # Extract timestep from the filename
        solid_file = vtkFile.replace('case_', 'solids_').replace('VTK/', 'VTK/solids/')
        
        if os.path.exists(solid_file):
            print(f"Loading solid structure from: {solid_file}")
            try:
                solid_reader = vtk.vtkUnstructuredGridReader()
                solid_reader.SetFileName(solid_file)
                solid_reader.Update()
                solid_data = solid_reader.GetOutput()
                
                # Try different approaches to get points
                if solid_data.GetNumberOfPoints() > 0:
                    # If the file has points directly
                    print(f"  Found {solid_data.GetNumberOfPoints()} points in solid structure")
                    solid_points = VN.vtk_to_numpy(solid_data.GetPoints().GetData())
                elif solid_data.GetNumberOfCells() > 0:
                    # If the file has cells, get cell centers
                    print(f"  Found {solid_data.GetNumberOfCells()} cells in solid structure")
                    solid_centers = vtk.vtkCellCenters()
                    solid_centers.SetInputData(solid_data)
                    solid_centers.Update()
                    centers_output = solid_centers.GetOutput()
                    if centers_output and centers_output.GetPoints():
                        solid_points = VN.vtk_to_numpy(centers_output.GetPoints().GetData())
                    else:
                        print("  Warning: Could not extract points from solid cell centers")
                else:
                    print("  Warning: Solid structure file contains no points or cells")
            except Exception as e:
                print(f"  Error processing solid structure: {e}")
        else:
            print(f"Warning: Solid structure file not found: {solid_file}")
    
    return {
        'points': fluid_points,
        'U': U,
        'Umag': Umag,
        'p': p,
        'Vcells': Vcells,
        'solid_points': solid_points
    }

def visualize_3d_velocity_field(data, time_step, plane='xy', z_slice=0, save=True):
    """Create velocity field visualization for a specific plane from 3D data
    Args:
        data: Dictionary containing the simulation data
        time_step: Time step of the simulation
        plane: Which plane to visualize ('xy', 'xz', or 'yz')
        z_slice: Position to slice the data (for xy-plane, this is the z-coordinate)
        save: Whether to save the plot
    """
    fig = plt.figure(figsize=(15, 6))
    
    # Get points and velocities
    points = data['points']
    U = data['U']
    Umag = data['Umag']
    solid_points = data.get('solid_points', None)
    
    # Select points near the chosen plane
    if plane == 'xy':
        mask = np.abs(points[:, 2] - z_slice) < 1e-5
        x, y = points[mask, 0], points[mask, 1]
        u, v = U[mask, 0], U[mask, 1]
        umag = Umag[mask]
        xlabel, ylabel = 'x [m]', 'y [m]'
        
        # Filter solid points if available
        if solid_points is not None:
            solid_mask = np.abs(solid_points[:, 2] - z_slice) < 1e-5
            solid_x, solid_y = solid_points[solid_mask, 0], solid_points[solid_mask, 1]
        else:
            solid_x, solid_y = None, None
            
    elif plane == 'xz':
        mask = np.abs(points[:, 1] - z_slice) < 1e-5
        x, y = points[mask, 0], points[mask, 2]
        u, v = U[mask, 0], U[mask, 2]
        umag = Umag[mask]
        xlabel, ylabel = 'x [m]', 'z [m]'
        
        # Filter solid points if available
        if solid_points is not None:
            solid_mask = np.abs(solid_points[:, 1] - z_slice) < 1e-5
            solid_x, solid_y = solid_points[solid_mask, 0], solid_points[solid_mask, 2]
        else:
            solid_x, solid_y = None, None
            
    else:  # yz plane
        mask = np.abs(points[:, 0] - z_slice) < 1e-5
        x, y = points[mask, 1], points[mask, 2]
        u, v = U[mask, 1], U[mask, 2]
        umag = Umag[mask]
        xlabel, ylabel = 'y [m]', 'z [m]'
        
        # Filter solid points if available
        if solid_points is not None:
            solid_mask = np.abs(solid_points[:, 0] - z_slice) < 1e-5
            solid_x, solid_y = solid_points[solid_mask, 1], solid_points[solid_mask, 2]
        else:
            solid_x, solid_y = None, None
    
    # Plot 1: Background scatter showing porous structure
    plt.subplot(121)
    scatter = plt.scatter(x, y, c=umag, s=2, cmap='viridis', alpha=0.8)
    
    # Add solid points if available
    if solid_points is not None and solid_x is not None and len(solid_x) > 0:
        plt.scatter(solid_x, solid_y, c='gray', s=1, alpha=0.6, label='Solid')
    
    plt.title(f'Time {time_step} - Velocity Magnitude Distribution\n{plane.upper()} Plane')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis('equal')
    plt.colorbar(scatter, label='Velocity Magnitude [m/s]')
    
    # Plot 2: Vector plot with contour background
    plt.subplot(122)
    # Use tricontourf for continuous velocity field visualization
    from matplotlib.tri import Triangulation
    triang = Triangulation(x, y)
    contour = plt.tricontourf(triang, umag, levels=20, cmap='viridis', alpha=0.9)
    
    # Add velocity vectors
    downsample = max(1, len(x) // 1000)
    plt.quiver(x[::downsample], y[::downsample],
              u[::downsample], v[::downsample],
              color='white', scale=15, width=0.002)
    
    # Add solid points if available
    if solid_points is not None and solid_x is not None and len(solid_x) > 0:
        plt.scatter(solid_x, solid_y, c='gray', s=1, alpha=0.6, label='Solid')
    
    plt.title(f'Time {time_step} - Flow Field with Vectors\n{plane.upper()} Plane')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis('equal')
    plt.colorbar(contour, label='Velocity Magnitude [m/s]')
    
    plt.tight_layout()
    if save:
        plt.savefig(f'figures_3d/velocity_field_{plane}_t{time_step}.png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def visualize_3d_pressure_field(data, time_step, save=True):
    """Create 3D pressure field visualization"""
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot colored by pressure
    scatter = ax.scatter(data['points'][:, 0],
                        data['points'][:, 1],
                        data['points'][:, 2],
                        c=data['p'],
                        cmap='viridis',
                        s=10)
    
    plt.colorbar(scatter, label='Pressure [Pa]')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    ax.set_title(f'3D Pressure Field - Time {time_step}')
    
    if save:
        plt.savefig(f'figures_3d/pressure_field_3d_t{time_step}.png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

#%%
data_0 = process_vtk_3d('VTK/case_0.vtk', 0, include_solids=True)
data_1000 = process_vtk_3d('VTK/case_1000.vtk', 1000, include_solids=True)

#%%
visualize_3d_velocity_field(data_0, 0, plane='xy', z_slice=0)

# #%%
# def visualize_3d_streamlines(data, time_step, save=True):
#     """Create 3D streamlines visualization using VTK"""
#     # Create a grid for streamline seeding
#     bounds = [
#         np.min(data['points'][:, 0]), np.max(data['points'][:, 0]),
#         np.min(data['points'][:, 1]), np.max(data['points'][:, 1]),
#         np.min(data['points'][:, 2]), np.max(data['points'][:, 2])
#     ]
    
#     # Create VTK points and vectors
#     points_vtk = vtk.vtkPoints()
#     vectors_vtk = vtk.vtkFloatArray()
#     vectors_vtk.SetNumberOfComponents(3)
    
#     for i in range(len(data['points'])):
#         points_vtk.InsertNextPoint(data['points'][i])
#         vectors_vtk.InsertNextTuple3(data['U'][i, 0], data['U'][i, 1], data['U'][i, 2])
    
#     # Create unstructured grid
#     grid = vtk.vtkUnstructuredGrid()
#     grid.SetPoints(points_vtk)
#     grid.GetPointData().SetVectors(vectors_vtk)
    
#     # Create streamlines
#     streamline = vtk.vtkStreamTracer()
#     streamline.SetInputData(grid)
    
#     # Set seed points
#     seeds = vtk.vtkPointSource()
#     seeds.SetNumberOfPoints(100)
#     seeds.SetCenter((bounds[0] + bounds[1])/2, bounds[2], (bounds[4] + bounds[5])/2)
#     seeds.SetRadius(min(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])/4)
#     seeds.Update()
    
#     streamline.SetSourceConnection(seeds.GetOutputPort())
#     streamline.SetMaximumPropagation(200)
#     streamline.SetInitialIntegrationStep(0.1)
#     streamline.SetIntegrationDirectionToBoth()
#     streamline.Update()
    
#     # Create streamline actor
#     mapper = vtk.vtkPolyDataMapper()
#     mapper.SetInputConnection(streamline.GetOutputPort())
    
#     actor = vtk.vtkActor()
#     actor.SetMapper(mapper)
#     actor.GetProperty().SetColor(1, 1, 1)
    
#     # Create renderer and window
#     renderer = vtk.vtkRenderer()
#     renderer.AddActor(actor)
#     renderer.SetBackground(0.1, 0.1, 0.1)
    
#     window = vtk.vtkRenderWindow()
#     window.AddRenderer(renderer)
#     window.SetSize(800, 800)
    
#     # Save to image
#     if save:
#         w2i = vtk.vtkWindowToImageFilter()
#         w2i.SetInput(window)
#         w2i.Update()
        
#         writer = vtk.vtkPNGWriter()
#         writer.SetFileName(f'figures_3d/streamlines_3d_t{time_step}.png')
#         writer.SetInputConnection(w2i.GetOutputPort())
#         writer.Write()
#     else:
#         window.Render()

# #%%
# def create_3d_visualizations():
#     # Process both time steps
#     data_0 = process_vtk_3d('VTK/case_0.vtk', 0)
#     data_537 = process_vtk_3d('VTK/case_537.vtk', 537)
    
#     # Create visualizations for both time steps
#     for time_step, data in [(0, data_0), (537, data_537)]:
#         print(f"\nCreating visualizations for time step {time_step}...")
        
#         # Velocity field visualization
#         visualize_3d_velocity_field(data, time_step)
        
#         # Pressure field visualization
#         visualize_3d_pressure_field(data, time_step)
        
#         # Streamlines visualization
#         visualize_3d_streamlines(data, time_step)
    
#     print("\nAll 3D visualizations have been saved in the 'figures_3d' directory.")

# if __name__ == "__main__":
#     create_3d_visualizations()

# %%
