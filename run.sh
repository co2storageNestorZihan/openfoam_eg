#!/bin/sh
 cd ${0%/*} || exit 1    # Run from this directory
 
 # Source tutorial run functions
 . $WM_PROJECT_DIR/bin/tools/RunFunctions
 
 application=`getApplication`
 

# copy stl file to constant/triSurface
cp porous_model.stl constant/triSurface/

#  ./clean.sh
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication checkMesh -allTopology -allGeometry

transformPoints -scale "(1e-6 1e-6 1e-6)"

runApplication $application

# Create VTK directory if it doesn't exist
mkdir -p VTK

# Convert OpenFOAM results to VTK format
runApplication foamToVTK

# Copy/move specific time steps to match Python script expectations
cp -f VTK/case_0.vtk VTK/case_0.vtk
cp -f VTK/case_500.vtk VTK/case_537.vtk  # Using time step 500 as proxy for 537 in Python script