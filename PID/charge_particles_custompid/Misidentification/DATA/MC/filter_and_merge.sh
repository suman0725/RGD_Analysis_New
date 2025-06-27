#!/bin/bash

# Base directory to search for all r_ttest3S.hipo files
BASE_DIR="/lustre24/expphy/volatile/clas12/dmat/test/febCarbonboth"

# Output directory for filtered files
OUT_DIR="/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/DATA/MC/temp_filtered"
MERGED_FILE="/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/DATA/MC/out.hipo"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Find all r_ttest3S.hipo files and process them
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
        # Run the filtering command
        hipo-utils -filter -b "MC::Particle","REC::Particle","REC::Cherenkov","REC::Scintillator","REC::Calorimeter" -o "$out_file" "$filepath"
    fi
done

# Merge all filtered files into one
echo "Merging filtered files into $MERGED_FILE"
hipo-utils -merge "$OUT_DIR"/filtered_*.hipo -o "$MERGED_FILE"

echo "Done!"
