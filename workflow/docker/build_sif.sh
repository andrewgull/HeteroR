#!/usr/bin/bash

CONTAINERS=( "assembly" "annotation" "biopython" "biostrings" "card_rgi" "default" "rscripts" )
OUTPUT_DIR="resources/apptainer"

# Check if the output directory exists, if not, create it
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Directory '$OUTPUT_DIR' does not exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi

for C in "${CONTAINERS[@]}"; do
    echo "Building SIF for container: $C"
    
    # Build Docker image
    docker build -t ${C} workflow/docker/${C}
    
    # Save Docker image to tar archive
    docker save -o workflow/docker/${C}.tar ${C}
    
    # Build Apptainer/Singularity image from Docker archive
    singularity build ${OUTPUT_DIR}/${C}.sif docker-archive://workflow/docker/${C}.tar
    
    # Clean up tar file after successful build
    rm -f workflow/docker/${C}.tar
    
    echo "Successfully built ${OUTPUT_DIR}/${C}.sif"
done

echo "All containers built successfully."