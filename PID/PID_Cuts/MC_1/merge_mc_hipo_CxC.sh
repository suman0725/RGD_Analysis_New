#!/bin/bash

# Base directory to search for all r_ttest3S.hipo files
BASE_DIR="/lustre24/expphy/volatile/clas12/dmat/clean/CCfullob"

# Output directory for filtered files
OUT_DIR="/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/temp_filtered"
MERGED_FILE="/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/out.hipo"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Find and process all r_ttest3S.hipo files
find "$BASE_DIR" -type f -name "r_ttest3S.hipo" | while read -r filepath; do 
    # Extract parent and grandparent directories for unique output name
    parent_dir=$(basename "$(dirname "$filepath")")
    grandparent_dir=$(basename "$(dirname "$(dirname "$filepath")")")

    # Unique output file
    out_file="$OUT_DIR/filtered_${grandparent_dir}_${parent_dir}.hipo"

    echo "Filtering $filepath -> $out_file"
    
    if [[ -f "$out_file" ]]; then
        echo "Skipping: $out_file already exists."
    else
        # Run the filtering command with all required banks
        hipo-utils -filter -b "MC::Particle,MC::Header,MC::Event,MC::Lund,RUN::config,REC::Particle,REC::Scintillator,REC::Cherenkov,REC::Event,REC::Calorimeter" -o "$out_file" "$filepath"
    fi
done

# Merge all filtered files into one
if ls "$OUT_DIR"/filtered_*.hipo >/dev/null 2>&1; then
    echo "Merging filtered files into $MERGED_FILE"
    hipo-utils -merge "$OUT_DIR"/filtered_*.hipo -o "$MERGED_FILE"
else
    echo "No filtered files found in $OUT_DIR. Skipping merge."
    exit 1
fi

echo "Done!"