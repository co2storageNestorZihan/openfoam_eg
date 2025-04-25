#!/bin/bash

# Clean script for OpenFOAM simulation directory

echo "Cleaning OpenFOAM simulation files..."

# Handle the 0 directory - keep only necessary files
if [ -d "0" ]; then
    echo "Cleaning 0 directory..."
    # Looking at the current structure, we need to preserve p 
    # Don't delete U if it exists as it was flagged as needed
    cd 0 2>/dev/null
    for file in *; do
        # Keep p and U, remove others like cellLevel, pointLevel, V
        if [[ "$file" != "p" && "$file" != "U" ]]; then
            echo "Removing 0/$file"
            rm -f "$file"
        fi
    done
    cd ..
fi

# Remove other time directories
for dir in $(find . -maxdepth 1 -type d -name "[0-9]*" | grep -v "^./0$"); do
    echo "Removing time directory: $dir"
    rm -rf "$dir"
done

# Keep only necessary files in the system directory
cd system 2>/dev/null
for file in *; do
    if [[ "$file" != "controlDict" && "$file" != "fvSchemes" && "$file" != "fvSolution" && \
          "$file" != "snappyHexMeshDict" && "$file" != "blockMeshDict" && "$file" != "meshQualityDict" ]]; then
        echo "Removing system/$file"
        rm -f "$file"
    fi
done
cd ..

# Keep only necessary files in the constant directory
cd constant 2>/dev/null
for file in *; do
    if [[ "$file" != "transportProperties" && "$file" != "turbulenceProperties" ]]; then
        echo "Removing constant/$file"
        rm -rf "$file"
    fi
done
cd ..

# Clean up processor directories
rm -rf processor*

# Remove OpenFOAM log files
rm -f log.*
rm -f *log
rm -f *.log

# Handle STL files - keep only those in geometry directory
echo "Cleaning STL files except in geometry directory..."
find . -path "./geometry" -prune -o -name "*.stl" -exec rm -f {} \;

# Remove ParaView files
rm -f *.foam
rm -f *.OpenFOAM

# Keep essential files in geometry directory
if [ -d "geometry" ]; then
    cd geometry 2>/dev/null
    for file in *; do
        if [[ "$file" != "porousModel.png" ]]; then
            # Keep STL files that might be input geometry
            if [[ "$file" != *.stl ]]; then
                echo "Removing geometry/$file"
                rm -f "$file"
            fi
        fi
    done
    cd ..
fi

# Preserve all Python, YAML, MD, and TXT files as per instructions
echo "Preserving Python, YAML, MD, and TXT files..."

echo "Clean completed."
