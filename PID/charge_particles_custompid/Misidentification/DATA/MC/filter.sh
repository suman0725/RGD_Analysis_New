#!/bin/bash

# Purpose: Filter MC::Particle from 5 r_ttest3S.hipo files and merge
BASE_DIR="/lustre24/expphy/volatile/clas12/dmat/test/febCarbonboth"
OUTPUT_DIR="/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/DATA/MC"
OUTPUT_FILE="$OUTPUT_DIR/filtered_r_ttest3S_all.hipo"
TEMP_DIR="$OUTPUT_DIR/temp"
LOG_FILE="$OUTPUT_DIR/filter_hipo.log"
HIPO_UTILS="/u/scigroup/cvmfs/hallb/clas12/sw/noarch/coatjava/11.1.1/bin/hipo-utils"

# Create directories and clear log
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"
> "$LOG_FILE"

# Collect 5 HIPO files
HIPO_FILES=()
COUNT=0
for SUBDIR in "$BASE_DIR"/*; do
    HIPO_FILE="$SUBDIR/sidis_mc-master/r_ttest3S.hipo"
    if [ -f "$HIPO_FILE" ]; then
        "$HIPO_UTILS" -test "$HIPO_FILE" >> "$LOG_FILE" 2>&1
        if [ $? -eq 0 ]; then
            HIPO_FILES+=("$HIPO_FILE")
            ((COUNT++))
            echo "Readable: $HIPO_FILE" >> "$LOG_FILE"
        fi
        if [ $COUNT -ge 5 ]; then
            break
        fi
    fi
done

if [ ${#HIPO_FILES[@]} -eq 0 ]; then
    echo "No readable files found" >> "$LOG_FILE"
    exit 1
fi

echo "Found ${#HIPO_FILES[@]} readable files" >> "$LOG_FILE"

# Filter MC::Particle
VALID_TEMP_FILES=()
for i in "${!HIPO_FILES[@]}"; do
    HIPO_FILE="${HIPO_FILES[$i]}"
    TEMP_OUTPUT="$TEMP_DIR/filtered_$i.hipo"
    echo "Filtering $HIPO_FILE to $TEMP_OUTPUT" >> "$LOG_FILE"
    "$HIPO_UTILS" -filter -b MC::Particle -s true -o "$TEMP_OUTPUT" "$HIPO_FILE" >> "$LOG_FILE" 2>&1
    if [ $? -eq 0 ] && [ -f "$TEMP_OUTPUT" ]; then
        echo "Success: $TEMP_OUTPUT created" >> "$LOG_FILE"
        VALID_TEMP_FILES+=("$TEMP_OUTPUT")
    else
        echo "Failed to process $HIPO_FILE" >> "$LOG_FILE"
    fi
done

# Merge
if [ ${#VALID_TEMP_FILES[@]} -gt 0 ]; then
    echo "Merging ${#VALID_TEMP_FILES[@]} files to $OUTPUT_FILE" >> "$LOG_FILE"
    "$HIPO_UTILS" -merge -o "$OUTPUT_FILE" "${VALID_TEMP_FILES[@]}" >> "$LOG_FILE" 2>&1
    if [ $? -eq 0 ]; then
        echo "Output file created: $OUTPUT_FILE" >> "$LOG_FILE"
    else
        echo "Error merging files" >> "$LOG_FILE"
        exit 1
    fi
else
    echo "No valid files to merge" >> "$LOG_FILE"
    exit 1
fi

# Clean up
rm -rf "$TEMP_DIR"
echo "Done" >> "$LOG_FILE"
