import vtk
import numpy as np
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt
import os
from typing import List, Dict, Optional, Tuple, Any

# --- Configuration Constants ---
# These can be overridden when calling functions if needed
DEFAULT_VTK_DIR = 'VTK'
DEFAULT_CT_DATA_FILE = 'large3d/subvolume.npy'
DEFAULT_VOXEL_RESOLUTION = 3.0035e-6 # meters/voxel (MUST MATCH SIMULATION SCALE)
# --- End Configuration Constants ---


def load_ct_data(ct_file_path: str = DEFAULT_CT_DATA_FILE) -> Optional[Tuple[np.ndarray, Tuple[int, int, int]]]:
    """Loads the 3D CT data from a .npy file."""
    print(f"Loading CT data from {ct_file_path}...")
    if not os.path.exists(ct_file_path):
        print(f"  Error: CT data file not found.")
        return None
    try:
        ct_data = np.load(ct_file_path)
        print(f"  CT data shape: {ct_data.shape}")
        return ct_data, ct_data.shape
    except Exception as e:
        print(f"  Error loading CT data: {e}")
        return None

def get_ct_slice(
    ct_data: np.ndarray,
    slice_dim: str,
    slice_frac_or_index: float | int,
    voxel_resolution: float = DEFAULT_VOXEL_RESOLUTION
) -> Optional[Dict[str, Any]]:
    """
    Extracts a 2D slice from the 3D CT data.

    Args:
        ct_data: The loaded 3D numpy array (CT data).
        slice_dim: Dimension to slice ('X', 'Y', or 'Z').
        slice_frac_or_index: Fractional position (0.0-1.0) or integer index for the slice.
        voxel_resolution: Physical size of one voxel in meters.

    Returns:
        A dictionary containing slice information, or None on error.
        Keys: 'ct_slice_2d', 'slice_coord', 'slice_axis_index', 'xlabel',
              'ylabel', 'plot_extent', 'slice_dim_name', 'slice_index'
    """
    Nx_ct, Ny_ct, Nz_ct = ct_data.shape
    Lx = Nx_ct * voxel_resolution
    Ly = Ny_ct * voxel_resolution
    Lz = Nz_ct * voxel_resolution

    slice_dim = slice_dim.upper()
    max_indices = {'X': Nx_ct - 1, 'Y': Ny_ct - 1, 'Z': Nz_ct - 1}

    if slice_dim not in max_indices:
        print(f"Error: Invalid slice_dim '{slice_dim}'. Choose 'X', 'Y', or 'Z'.")
        return None

    # Determine slice index and coordinate
    if isinstance(slice_frac_or_index, float):
        if not (0.0 <= slice_frac_or_index <= 1.0):
            print("Error: slice_frac must be between 0.0 and 1.0.")
            return None
        slice_index = int(max_indices[slice_dim] * slice_frac_or_index)
    elif isinstance(slice_frac_or_index, int):
        slice_index = slice_frac_or_index
    else:
         print("Error: slice_frac_or_index must be a float (0-1) or an int.")
         return None

    if not (0 <= slice_index <= max_indices[slice_dim]):
        print(f"Error: slice_index {slice_index} is out of bounds for dimension {slice_dim} (0-{max_indices[slice_dim]}).")
        return None

    slice_coord = (slice_index + 0.5) * voxel_resolution

    # Extract slice and define plotting parameters
    if slice_dim == 'Z':
        ct_slice_2d = ct_data[:, :, slice_index].T # Transpose for imshow
        slice_axis_index = 2
        xlabel, ylabel = 'x [m]', 'y [m]'
        plot_extent = [0, Lx, 0, Ly]
    elif slice_dim == 'Y':
        ct_slice_2d = ct_data[:, slice_index, :].T
        slice_axis_index = 1
        xlabel, ylabel = 'x [m]', 'z [m]'
        plot_extent = [0, Lx, 0, Lz]
    elif slice_dim == 'X':
        ct_slice_2d = ct_data[slice_index, :, :].T
        slice_axis_index = 0
        xlabel, ylabel = 'y [m]', 'z [m]'
        plot_extent = [0, Ly, 0, Lz]

    print(f"Selected CT slice: {slice_dim} = {slice_index} (physical coord: {slice_coord:.3e} m)")

    return {
        'ct_slice_2d': ct_slice_2d,
        'slice_coord': slice_coord,
        'slice_axis_index': slice_axis_index,
        'xlabel': xlabel,
        'ylabel': ylabel,
        'plot_extent': plot_extent,
        'slice_dim_name': slice_dim,
        'slice_index': slice_index
    }


def load_sim_data(vtk_file_path: str) -> Optional[Dict[str, np.ndarray]]:
    """
    Loads simulation data (coordinates, U, p) from a single VTK file.

    Returns: Dictionary with keys 'x', 'y', 'z', 'Umag', 'p', 'Ux', 'Uy', 'Uz',
             or None on failure.
    """
    print(f"Processing {vtk_file_path}...")
    if not os.path.exists(vtk_file_path):
        print(f"  Error: VTK file not found.")
        return None

    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file_path)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()

    if data is None or data.GetNumberOfCells() == 0:
        print(f"  Error: Failed to read data or no cells found.")
        return None

    # Check for required arrays
    if data.GetCellData().GetArray('U') is None or data.GetCellData().GetArray('p') is None:
        print(f"  Error: 'U' or 'p' array not found in cell data.")
        return None

    # Extract cell centers
    try:
        cell_centers_filter = vtk.vtkCellCenters()
        cell_centers_filter.SetInputData(data)
        cell_centers_filter.Update()
        centers_data = cell_centers_filter.GetOutput()
        points_vtk = centers_data.GetPoints().GetData()
        coords = VN.vtk_to_numpy(points_vtk)
        x_coords, y_coords, z_coords = coords[:,0], coords[:,1], coords[:,2]

        # Extract velocity and pressure
        U = VN.vtk_to_numpy(data.GetCellData().GetArray('U'))
        p = VN.vtk_to_numpy(data.GetCellData().GetArray('p'))
        Ux, Uy, Uz = U[:,0], U[:,1], U[:,2]
        Umag = np.sqrt(Ux**2 + Uy**2 + Uz**2)

        print(f"  Successfully processed {data.GetNumberOfCells()} cells.")
        return {'x': x_coords, 'y': y_coords, 'z': z_coords,
                'Umag': Umag, 'p': p, 'Ux': Ux, 'Uy': Uy, 'Uz': Uz}
    except Exception as e:
        print(f"  Error extracting data from VTK: {e}")
        return None


def plot_slice_comparison(
    ct_slice_info: Dict[str, Any],
    sim_data: Dict[str, np.ndarray],
    time_step: int,
    field_to_plot: str,
    voxel_resolution: float = DEFAULT_VOXEL_RESOLUTION,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    s: int = 2, # Scatter point size
    cmap: str = 'viridis'
) -> Optional[plt.Figure]:
    """
    Generates a 1x2 plot comparing a CT slice with simulation data for a specific field.

    Args:
        ct_slice_info: Dictionary returned by get_ct_slice.
        sim_data: Dictionary returned by load_sim_data for the specific time_step.
        time_step: The simulation time step being plotted.
        field_to_plot: The simulation field to plot ('Umag' or 'p').
        voxel_resolution: Physical size of one voxel in meters.
        vmin, vmax: Optional min/max values for the simulation data color scale.
        s: Size of scatter points for simulation data.
        cmap: Colormap for simulation data.

    Returns:
        A matplotlib Figure object, or None if plotting fails.
    """
    if field_to_plot not in ['Umag', 'p']:
        print("Error: field_to_plot must be 'Umag' or 'p'.")
        return None
    if field_to_plot not in sim_data:
         print(f"Error: Field '{field_to_plot}' not found in simulation data.")
         return None

    slice_coord = ct_slice_info['slice_coord']
    slice_axis_index = ct_slice_info['slice_axis_index']
    slice_tolerance = voxel_resolution / 2.0
    slice_dim_name = ct_slice_info['slice_dim_name']
    slice_index = ct_slice_info['slice_index']

    # Filter simulation data near the selected slice
    coords_all = np.stack((sim_data['x'], sim_data['y'], sim_data['z']), axis=-1)
    slice_mask = np.abs(coords_all[:, slice_axis_index] - slice_coord) < slice_tolerance

    if not np.any(slice_mask):
        print(f"Warning: No simulation cells found near slice at {slice_dim_name}={slice_coord:.3e} m for time {time_step}.")
        return None

    # Determine plotting coordinates based on slice dimension
    if slice_dim_name == 'Z':
        x_plot = sim_data['x'][slice_mask]
        y_plot = sim_data['y'][slice_mask]
    elif slice_dim_name == 'Y':
        x_plot = sim_data['x'][slice_mask]
        y_plot = sim_data['z'][slice_mask]
    elif slice_dim_name == 'X':
        x_plot = sim_data['y'][slice_mask]
        y_plot = sim_data['z'][slice_mask]

    sim_field_slice = sim_data[field_to_plot][slice_mask]

    # Create figure
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5.5))
    ax_ct, ax_sim = axes[0], axes[1]

    # Plot CT Slice
    ax_ct.imshow(ct_slice_info['ct_slice_2d'], cmap='gray', origin='lower',
                 extent=ct_slice_info['plot_extent'])
    ax_ct.set_title(f'CT Slice ({slice_dim_name}={slice_index})')
    ax_ct.set_xlabel(ct_slice_info['xlabel'])
    ax_ct.set_ylabel(ct_slice_info['ylabel'])
    ax_ct.set_aspect('equal')

    # Plot Simulation Field Slice
    if field_to_plot == 'Umag' and vmin is None:
        vmin = 0 # Sensible default for velocity magnitude

    sc = ax_sim.scatter(x_plot, y_plot, c=sim_field_slice, s=s, cmap=cmap, vmin=vmin, vmax=vmax)
    field_label = 'Velocity Magnitude [m/s]' if field_to_plot == 'Umag' else 'Pressure [Pa]'
    ax_sim.set_title(f'{field_label.split(" ")[0]} Field (Time = {time_step})')
    ax_sim.set_xlabel(ct_slice_info['xlabel'])
    ax_sim.set_ylabel(ct_slice_info['ylabel'])
    ax_sim.set_xlim(ct_slice_info['plot_extent'][0], ct_slice_info['plot_extent'][1])
    ax_sim.set_ylim(ct_slice_info['plot_extent'][2], ct_slice_info['plot_extent'][3])
    ax_sim.set_aspect('equal')
    fig.colorbar(sc, ax=ax_sim, label=field_label)

    plt.tight_layout()
    return fig


def plot_velocity_histograms(
    time_steps: List[int],
    vtk_dir: str = DEFAULT_VTK_DIR
) -> Optional[plt.Figure]:
    """
    Generates a histogram plot comparing velocity magnitude distributions
    for multiple time steps.

    Args:
        time_steps: A list of integer time steps to compare.
        vtk_dir: The directory containing the VTK files (e.g., 'VTK').

    Returns:
        A matplotlib Figure object, or None if no data could be plotted.
    """
    print("\nGenerating Velocity Histograms...")
    fig, ax = plt.subplots(figsize=(8, 6))
    has_hist_data = False

    all_umag_data = [] # Collect all data to determine overall range if needed

    # Load data first
    sim_data_all_times = {}
    for t in time_steps:
        vtk_file = os.path.join(vtk_dir, f'case_{t}.vtk')
        data = load_sim_data(vtk_file)
        if data is not None and 'Umag' in data:
            sim_data_all_times[t] = data['Umag']
            all_umag_data.append(data['Umag'])
        else:
            print(f"Warning: Could not load Umag data for time step {t}. Skipping histogram.")

    if not sim_data_all_times:
        print("Error: No valid velocity magnitude data found for any specified time step.")
        plt.close(fig)
        return None

    # Determine common range or plot individually
    # For simplicity here, we plot individually scaled histograms
    for t, umag_data in sim_data_all_times.items():
        if len(umag_data) > 0:
             ax.hist(umag_data, bins=50, alpha=0.7, label=f'Time {t}', density=True)
             has_hist_data = True

    if not has_hist_data:
         print("Error: Although files were loaded, no valid data points found for histograms.")
         plt.close(fig)
         return None

    ax.set_title('Velocity Magnitude Distribution (Full Domain)')
    ax.set_xlabel('Velocity Magnitude [m/s]')
    ax.set_ylabel('Density')
    ax.legend()
    # Set y-axis to log scale only if data warrants it (avoid errors if all counts are zero)
    y_lims = ax.get_ylim()
    if y_lims[1] > y_lims[0] > 0: # Check if range is valid and positive
         ax.set_yscale('log')
    else:
        print("Note: Could not set y-axis to log scale (likely due to data range).")

    plt.tight_layout()
    return fig

# --- Example Usage (for testing in a standard Python script) ---
# In a Jupyter notebook, you would typically call the functions individually.
if __name__ == "__main__":

    OUTPUT_DIR_MAIN = 'figures_vis_interactive_test'
    if not os.path.exists(OUTPUT_DIR_MAIN):
        os.makedirs(OUTPUT_DIR_MAIN)

    # 1. Load CT
    ct_result = load_ct_data()
    if ct_result is None:
        exit()
    ct_data_main, ct_shape_main = ct_result

    # 2. Define Slice Parameters
    slice_dimension_main = 'Z'
    slice_position_main = 0.5 # Middle slice fraction

    # 3. Get CT Slice Info
    ct_slice_info_main = get_ct_slice(ct_data_main, slice_dimension_main, slice_position_main)
    if ct_slice_info_main is None:
        exit()

    # 4. Define Time Steps and Load Sim Data
    time_steps_main = [0, 500] # Use available time steps
    sim_data_cache = {}
    for t in time_steps_main:
         vtk_file = os.path.join(DEFAULT_VTK_DIR, f'case_{t}.vtk')
         sim_data_cache[t] = load_sim_data(vtk_file)

    # 5. Plot Velocity Comparison for each time step
    print("\n--- Generating Example Plots ---")
    for t in time_steps_main:
        if sim_data_cache.get(t): # Check if data was loaded successfully
            print(f"\nPlotting comparison for Time={t}, Field=Umag")
            fig_vel = plot_slice_comparison(ct_slice_info_main, sim_data_cache[t], t, 'Umag', cmap='viridis')
            if fig_vel:
                 fig_vel_path = os.path.join(OUTPUT_DIR_MAIN, f'example_CT_vs_Velocity_T{t}.png')
                 fig_vel.savefig(fig_vel_path, dpi=150)
                 plt.close(fig_vel) # Close figure to avoid displaying in script mode
                 print(f"Saved example plot: {fig_vel_path}")

            print(f"\nPlotting comparison for Time={t}, Field=p")
            # Determine pressure range for consistency (optional, but good practice)
            all_p = [sd['p'] for sd in sim_data_cache.values() if sd is not None]
            p_min_glob = np.min(np.concatenate(all_p)) if all_p else 0
            p_max_glob = np.max(np.concatenate(all_p)) if all_p else 1

            fig_p = plot_slice_comparison(ct_slice_info_main, sim_data_cache[t], t, 'p', cmap='plasma', vmin=p_min_glob, vmax=p_max_glob)
            if fig_p:
                 fig_p_path = os.path.join(OUTPUT_DIR_MAIN, f'example_CT_vs_Pressure_T{t}.png')
                 fig_p.savefig(fig_p_path, dpi=150)
                 plt.close(fig_p)
                 print(f"Saved example plot: {fig_p_path}")
        else:
            print(f"\nSkipping plots for Time={t} due to load error.")


    # 6. Plot Velocity Histograms
    fig_hist = plot_velocity_histograms(time_steps_main)
    if fig_hist:
        fig_hist_path = os.path.join(OUTPUT_DIR_MAIN, 'example_Velocity_Histograms.png')
        fig_hist.savefig(fig_hist_path, dpi=150)
        plt.close(fig_hist)
        print(f"\nSaved example histogram: {fig_hist_path}")

    print("\nExample usage finished.")
