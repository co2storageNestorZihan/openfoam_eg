#!/bin/bash
cd ${0%/*} || exit 1 # run from this directory

echo "Cleaning simulation outputs while preserving configuration files..."

# Remove time step directories (but not 0 directory which contains configuration)
for d in [0-9]* [0-9]*[0-9]; do
    # Skip the 0 directory which contains configuration files
    if [ "$d" != "0" ]; then
        if [ -d "$d" ]; then
            echo "Removing time step directory: $d"
            rm -rf "$d"
        fi
    fi
done

# Remove log files
echo "Removing log files..."
rm -f log.* 

# Clean VTK directory but keep the structure
if [ -d "VTK" ]; then
    echo "Cleaning VTK directory..."
    rm -f VTK/*.vtk
    # Keep directory structure
    for d in VTK/*; do
        if [ -d "$d" ]; then
            echo "Keeping directory structure: $d"
            rm -f "$d"/*
        fi
    done
fi

# Remove generated figure files
if [ -d "figures" ]; then
    echo "Removing generated figures..."
    rm -rf figures
fi
rm -f velocity_maps.png

# Clean postProcessing directory if it exists
if [ -d "postProcessing" ]; then
    echo "Removing postProcessing directory..."
    rm -rf postProcessing
fi

# Clean dynamicCode directory if it exists
if [ -d "dynamicCode" ]; then
    echo "Removing dynamicCode directory..."
    rm -rf dynamicCode
fi

# Clean processor* directories if they exist
rm -rf processor*

# Clean OpenFOAM-specific temporary files
rm -f *.OpenFOAM
rm -f *.foam

# Clean mesh backup
rm -rf constant/polyMesh/sets
rm -rf constant/polyMesh/*Zones*
rm -rf constant/polyMesh/refinementHistory
rm -rf constant/extendedFeatureEdgeMesh

# But keep the essential mesh files
if [ -d "constant/polyMesh" ]; then
    echo "Preserving essential mesh files in constant/polyMesh..."
fi

echo "Cleaning completed."