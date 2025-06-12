#!/bin/sh
 cd ${0%/*} || exit 1    # Run from this directory
 
 # Source tutorial run functions
 . $WM_PROJECT_DIR/bin/tools/RunFunctions
 
 application=`getApplication`
 

# copy stl file to constant/triSurface
cp porous_model.stl constant/triSurface/

# Read locationInMesh coordinates from external file
# Priority: 1) Command line argument, 2) Environment variable, 3) Default path
if [ "$1" != "" ]; then
    LOCATION_FILE="$1/locationInMesh.txt"
elif [ "$LOCATION_DATA_PATH" != "" ]; then
    LOCATION_FILE="$LOCATION_DATA_PATH/locationInMesh.txt"
else
    LOCATION_FILE="/Volumes/topaz_ssd4tb_macbook/porescale_data_small/bentheimer/idx_0/locationInMesh.txt"
fi

if [ -f "$LOCATION_FILE" ]; then
    echo "Reading locationInMesh coordinates from $LOCATION_FILE"
    
    # Read the three coordinates from the file
    X_COORD=$(sed -n '1p' "$LOCATION_FILE")
    Y_COORD=$(sed -n '2p' "$LOCATION_FILE")
    Z_COORD=$(sed -n '3p' "$LOCATION_FILE")
    
    echo "Using locationInMesh coordinates: ($X_COORD $Y_COORD $Z_COORD)"
    
    # Update the snappyHexMeshDict file with new coordinates
    sed -i.bak "s/locationInMesh.*$/locationInMesh ($X_COORD $Y_COORD $Z_COORD);/" system/snappyHexMeshDict
    
    echo "Updated snappyHexMeshDict with new locationInMesh coordinates"
else
    echo "Warning: $LOCATION_FILE not found. Using default locationInMesh coordinates."
fi

#  ./clean.sh
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication checkMesh -allTopology -allGeometry

transformPoints -scale "(3.0035e-6 3.0035e-6 3.0035e-6)"

runApplication $application

# Create VTK directory if it doesn't exist
mkdir -p VTK

# Convert OpenFOAM results to VTK format
runApplication foamToVTK

# Copy/move specific time steps to match Python script expectations
cp -f VTK/case_0.vtk VTK/case_0.vtk
cp -f VTK/case_500.vtk VTK/case_537.vtk  # Using time step 500 as proxy for 537 in Python script