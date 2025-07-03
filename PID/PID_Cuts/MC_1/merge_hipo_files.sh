#!/bin/bash

# Base directory to search for all r_ttest3S.hipo files
BASE_DIR="/lustre24/expphy/volatile/clas12/dmat/clean/CCfullob"

# Output directory for filtered files
OUT_DIR="/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/temp_filtered"
MERGED_FILE="/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/out.hipo"

# Option to clean up log files after successful run (set to 1 to enable cleanup)
CLEAN_LOGS=0

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Enable debugging
set -x

# Find and process all r_ttest3S.hipo files
find "$BASE_DIR" -type f -name "r_ttest3S.hipo" | while read -r filepath; do 
    # Extract parent and grandparent directories for unique output name
    parent_dir=$(basename "$(dirname "$filepath")")
    grandparent_dir=$(basename "$(dirname "$(dirname "$filepath")")")
    base_name=$(basename "$filepath" .hipo) # Strip .hipo extension

    # Unique output file
    out_file="$OUT_DIR/filtered_${grandparent_dir}_${parent_dir}_${base_name}.hipo"

    echo "Filtering $filepath -> $out_file"
    
    if [[ -f "$out_file" ]]; then
        echo "Skipping: $out_file already exists."
    else
        # Run the filtering command with all required banks and log errors
        hipo-utils -filter -b "MC::Particle,MC::Header,MC::Event,MC::Lund,RUN::config,REC::Particle,REC::Scintillator,REC::Cherenkov,REC::Sang,REC::Calorimeter" \
        -o "$out_file" "$filepath" 2> "$out_file.error.log"
        if [[ $? -ne 0 ]]; then
            echo "Error filtering $filepath. Check $out_file.error.log"
            mv "$out_file" "$out_file.corrupted" 2>/dev/null || true
            continue
        fi
        # Validate the filtered file (check file size > 100 KB as a heuristic)
        file_size=$(stat -c %s "$out_file" 2>/dev/null || echo 0)
        if [[ $file_size -lt 100000 ]]; then
            echo "Warning: $out_file is too small ($file_size bytes), likely corrupted. Skipping."
            mv "$out_file" "$out_file.corrupted" 2>/dev/null || true
            continue
        fi
    fi
done

# Merge all valid filtered files into one
if ls "$OUT_DIR"/filtered_*.hipo >/dev/null 2>&1; then
    echo "Merging filtered files into $MERGED_FILE"
    # Create a list of valid files to merge
    valid_files=()
    for file in "$OUT_DIR"/filtered_*.hipo; do
        # Verify file integrity before merging
        hipo-utils -info "$file" >/dev/null 2> "$file.info.log"
        if [[ $? -eq 0 ]]; then
            valid_files+=("$file")
        else
            echo "Warning: $file is invalid or corrupted. Excluding from merge. Check $file.info.log"
            mv "$file" "$file.corrupted" 2>/dev/null || true
        fi
    done
    if [[ ${#valid_files[@]} -eq 0 ]]; then
        echo "No valid filtered files found in $OUT_DIR. Skipping merge."
        exit 1
    fi
    echo "Merging ${#valid_files[@]} valid files"
    # Merge valid files
    hipo-utils -merge "${valid_files[@]}" -o "$MERGED_FILE" 2> "$MERGED_FILE.error.log"
    if [[ $? -ne 0 ]]; then
        echo "Error merging files. Check $MERGED_FILE.error.log"
        exit 1
    fi
else
    echo "No filtered files found in $OUT_DIR. Skipping merge."
    exit 1
fi

# Clean up log files if enabled and merge was successful
if [[ $CLEAN_LOGS -eq 1 ]]; then
    echo "Cleaning up log files"
    rm -f "$OUT_DIR"/*.error.log "$OUT_DIR"/*.info.log "$MERGED_FILE.error.log"
fi

echo "Done!"