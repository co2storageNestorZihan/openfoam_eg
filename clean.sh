#!/bin/bash
cd ${0%/*} || exit 1 # run from this directory

echo "Cleaning simulation outputs while preserving essential configuration files..."

# Files to keep (specify explicit patterns)
declare -a keep_patterns=(
    # 0 directory files
    "0/U"
    "0/V" 
    "0/p"
    # constant directory files
    "constant/transportProperties"
    "constant/turbulenceProperties"
    # geometry files
    "geometry/*.png"
    # system files
    "system/blockMeshDict"
    "system/controlDict"
    "system/fvSchemes"
    "system/fvSolution"
    "system/meshQualityDict"
    "system/snappyHexMeshDict"
    # scripts
    "clean.sh"
    "run.sh"
    # all Python files
    "*.py"
    # documentation and configuration files
    "*.md"
    "*.yaml"
    "*.yml"
    "*.txt"
)

# Create a temporary directory
temp_dir=$(mktemp -d)
echo "Using temporary directory: $temp_dir"

# Copy files to keep to temp directory, preserving directory structure
for pattern in "${keep_patterns[@]}"; do
    for file in $pattern; do
        if [ -f "$file" ]; then
            mkdir -p "$temp_dir/$(dirname "$file")"
            echo "Preserving: $file"
            cp "$file" "$temp_dir/$(dirname "$file")/"
        fi
    done
done

# Special handling for 0, constant, system, and geometry directories
# This ensures we only keep the directories themselves, not their contents
for dir in "0" "constant" "system" "geometry"; do
    if [ -d "$dir" ]; then
        mkdir -p "$temp_dir/$dir"
    fi
done

# Remove everything except .git directory
find . -mindepth 1 -not -path "./.git*" -not -path "./clean.sh" -not -path "$temp_dir*" -exec rm -rf {} \;

# Copy preserved files back
cp -r "$temp_dir"/* .

# Remove temporary directory
rm -rf "$temp_dir"

echo "Cleaning completed. Only essential files for simulation initialization are preserved."