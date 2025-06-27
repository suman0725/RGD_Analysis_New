#!/bin/bash


BASE_DIR="/lustre24/expphy/volatile/clas12/dmat/test/febCarbonboth"
OUTPUT_DIR="/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/DATA/MC"
OUTPUT_FILE="$OUTPUT_DIR/out.hipo"
TEMP_DIR="$OUTPUT_DIR/temp"
LOG_FILE="$OUTPUT_DIR/filter_hipo.log"
HIPO_UTILS="/u/scigroup/cvmfs/hallb/clas12/sw/noarch/clara/5.0.2_11.1.1/plugins/clas12/bin/hipo-utils"

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"
> "$LOG_FILE"

# Collect r_ttest3S.hipo files (limit to 5)
HIPO_FILES=($(find "$BASE_DIR" -type f -path "*/sidis_mc-master/r_ttest3S.hipo" | head -n 5))
if [ ${#HIPO_FILES[@]} -eq 0 ]; then
    echo "No r_ttest3S.hipo files found" >> "$LOG_FILE"
    exit 1
fi

echo "Found ${#HIPO_FILES[@]} files" >> "$LOG_FILE"

# Filter files
VALID_TEMP_FILES=()
for i in "${!HIPO_FILES[@]}"; do
    HIPO_FILE="${HIPO_FILES[$i]}"
    TEMP_OUTPUT="$TEMP_DIR/filtered_$i.hipo"
    echo "Filtering $HIPO_FILE" >> "$LOG_FILE"
    "$HIPO_UTILS" -filter -b "\"MC::Particle\",\"REC::Particle\",\"REC::Cherenkov\",\"REC::Scintillator\",\"REC::Calorimeter\"" -o "$TEMP_OUTPUT" "$HIPO_FILE" >> "$LOG_FILE" 2>&1
    if [ $? -eq 0 ] && [ -f "$TEMP_OUTPUT" ]; then
        VALID_TEMP_FILES+=("$TEMP_OUTPUT")
        echo "Success: $TEMP_OUTPUT created" >> "$LOG_FILE"
    else
        echo "Failed: $HIPO_FILE" >> "$LOG_FILE"
    fi
done

# Merge
if [ ${#VALID_TEMP_FILES[@]} -gt 0 ]; then
    echo "Merging ${#VALID_TEMP_FILES[@]} files to $OUTPUT_FILE" >> "$LOG_FILE"
    "$HIPO_UTILS" -merge -o "$OUTPUT_FILE" "${VALID_TEMP_FILES[@]}" >> "$LOG_FILE" 2>&1
    if [ $? -eq 0 ]; then
        echo "Output created: $OUTPUT_FILE" >> "$LOG_FILE"
    else
        echo "Merge failed" >> "$LOG_FILE"
        exit 1
    fi
else
    echo "No files to merge" >> "$LOG_FILE"
    exit 1
fi

rm -rf "$TEMP_DIR"
echo "Done" >> "$LOG_FILE"


