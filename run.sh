#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=$(getApplication)

# Create timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="simulation_log_${TIMESTAMP}.log"

# Function to log with timestamp
log_with_time() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to format seconds to human readable format
format_time() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))
    
    if [ $hours -gt 0 ]; then
        echo "${hours}h ${minutes}m ${seconds}s"
    elif [ $minutes -gt 0 ]; then
        echo "${minutes}m ${seconds}s"
    else
        echo "${seconds}s"
    fi
}

# Record overall start time
OVERALL_START=$(date +%s)
log_with_time "=== BATCH SIMULATION STARTED ==="

# Read the external folder absolute path from cfg.yaml
# Using grep and awk to parse the yaml file
EXT_DATA_DIR=$(grep 'root_data_dir:' cfg.yaml | awk '{print $2}')

if [ -z "$EXT_DATA_DIR" ] || [ ! -d "$EXT_DATA_DIR" ]; then
    log_with_time "ERROR: root_data_dir not found or is not a directory in cfg.yaml."
    exit 1
fi

log_with_time "External data directory: $EXT_DATA_DIR"

# Count total cases
TOTAL_CASES=$(find "$EXT_DATA_DIR" -maxdepth 1 -type d -name "idx_*" | wc -l)
log_with_time "Found $TOTAL_CASES cases to process"

CASE_COUNT=0
SUCCESSFUL_CASES=0

# Find all idx_* subfolders and loop through them
for CASE_DIR in "$EXT_DATA_DIR"/idx_*; do
    if [ -d "$CASE_DIR" ]; then
        CASE_COUNT=$((CASE_COUNT + 1))
        CASE_START=$(date +%s)
        
        echo "----------------------------------------------------"
        echo "Processing case $CASE_COUNT/$TOTAL_CASES: $(basename "$CASE_DIR")"
        echo "----------------------------------------------------"
        log_with_time "Started case $(basename "$CASE_DIR")"

        STL_FILE="$CASE_DIR/porous_model.stl"
        SNAPPY_DICT="$CASE_DIR/snappyHexMeshDict"

        if [ ! -f "$STL_FILE" ]; then
            echo "Warning: porous_model.stl not found in $CASE_DIR. Skipping..."
            continue
        fi

        if [ ! -f "$SNAPPY_DICT" ]; then
            echo "Warning: snappyHexMeshDict not found in $CASE_DIR. Skipping..."
            continue
        fi

        # Clean the case before running
        echo "Cleaning up previous run..."
        ./clean.sh
        echo "Done cleaning."

        # Copy stl file to the right folder
        echo "Copying $STL_FILE to constant/triSurface/"
        cp "$STL_FILE" constant/triSurface/

        # Copy snappyHexMeshDict to system (force overwrite)
        echo "Copying $SNAPPY_DICT to system/"
        cp -f "$SNAPPY_DICT" system/

        # Run simulation steps
        echo "Running blockMesh..."
        runApplication blockMesh
        echo "Running snappyHexMesh..."
        runApplication snappyHexMesh -overwrite
        echo "Running checkMesh..."
        runApplication checkMesh -allTopology -allGeometry

        echo "Running transformPoints..."
        transformPoints -scale "(3.0035e-6 3.0035e-6 3.0035e-6)"

        echo "Running solver: $application..."
        runApplication $application

        echo "Converting to VTK format..."
        runApplication foamToVTK

        # Copy specific VTK results to the idx folder
        RESULTS_DIR="$CASE_DIR/VTK_results"
        echo "Copying specific VTK results to $RESULTS_DIR"
        mkdir -p "$RESULTS_DIR"
        
        # Only copy case_0.vtk and case_500.vtk if they exist
        if [ -f "VTK/case_0.vtk" ]; then
            cp "VTK/case_0.vtk" "$RESULTS_DIR/"
            echo "Copied case_0.vtk"
        else
            echo "Warning: case_0.vtk not found"
        fi
        
        if [ -f "VTK/case_500.vtk" ]; then
            cp "VTK/case_500.vtk" "$RESULTS_DIR/"
            echo "Copied case_500.vtk"
        else
            echo "Warning: case_500.vtk not found"
        fi

        # Calculate case duration
        CASE_END=$(date +%s)
        CASE_DURATION=$((CASE_END - CASE_START))
        SUCCESSFUL_CASES=$((SUCCESSFUL_CASES + 1))
        
        echo "Case $(basename "$CASE_DIR") finished in $(format_time $CASE_DURATION). Results saved to $RESULTS_DIR"
        log_with_time "Case $(basename "$CASE_DIR") completed in $(format_time $CASE_DURATION)"
    fi
done

# Calculate total duration
OVERALL_END=$(date +%s)
TOTAL_DURATION=$((OVERALL_END - OVERALL_START))

echo "----------------------------------------------------"
echo "All simulations finished."
echo "Total runtime: $(format_time $TOTAL_DURATION)"
echo "Successful cases: $SUCCESSFUL_CASES/$CASE_COUNT"
echo "----------------------------------------------------"

log_with_time "=== BATCH SIMULATION COMPLETED ==="
log_with_time "Total runtime: $(format_time $TOTAL_DURATION)"
log_with_time "Successful cases: $SUCCESSFUL_CASES/$CASE_COUNT"
if [ $SUCCESSFUL_CASES -gt 0 ]; then
    log_with_time "Average time per case: $(format_time $((TOTAL_DURATION / SUCCESSFUL_CASES)))"
fi