# 2D Permeability Calculation Framework

This repository contains tools for calculating absolute permeability in 2D porous media using OpenFOAM simulation data.

Reproduction from: https://lruepke.github.io/HTF_lecture/summer2022/lectures/L03/FirstCase.html#making-the-mesh




## OpenFOAM Configuration

### Solver

The simulations use the `simpleFoam` solver from OpenFOAM, which:
- Implements the SIMPLE (Semi-Implicit Method for Pressure Linked Equations) algorithm
- Solves steady-state, incompressible Navier-Stokes equations
- Is suitable for laminar flows through porous media

### Boundary Conditions

The simulations are configured with:
- Fixed pressure boundary conditions at inlet (1 Pa) and outlet (0 Pa)
- No-slip conditions on solid walls
- Pressure gradient applied in the x-direction

## Data Pipeline

1. **Simulation**: OpenFOAM simulations produce VTK files containing velocity and pressure fields
2. **Data Processing**: The Python script processes the VTK files to extract velocity and pressure fields
3. **Permeability Calculation**: Darcy velocity is computed from cell volumes and velocities, then used to calculate permeability
4. **Visualization**: The script generates spatial distributions and histograms for analysis

## Permeability Calculation

The calculation follows Darcy's Law:
- Volume fluxes are calculated for each cell: qₓ = Ux × Vcell
- Total flux is calculated by summing individual cell fluxes
- Darcy velocity is obtained by dividing total flux by total domain volume (not just pore volume)
- Permeability is calculated using: kxx = μ·L·UDarcy / ΔP
- Results are converted to millidarcy (md) for reporting

## Prerequisites and Running the Pipeline

### Prerequisites

1. Install required Python packages using conda or miniforge:
```bash
conda env create -f environment.yaml
conda activate permeability-env
```

### Running the Simulation

1. First, run the image conversion script to generate VTI and STL files:
```bash
python img_convert.py
```

2. Run the simulation using the provided shell script:
```bash
./run.sh
```

3. Calculate permeability and generate visualizations:
```bash
python calculate_perm.py
```

4. Results will be printed to console and figures saved to the `figures/` directory
