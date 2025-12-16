#!/bin/bash

# nohup ./check_mergingstep.sh > nohup_check_mergingstep.out 2>&1 &

# [INFO] Script to check if the byte size of a merged file (CTGA)
# [INFO] is exactly equal to the sum of the byte sizes of its two components (CT + GA).
# [INFO] It iterates over specified sample directories and threshold levels (3, 4, 5).

# --- Configuration ---
# The base directory containing all sample directories (e.g., AID_KO_CIT, UM_TKO_CIT)
BASE_DIR="/30tb/home/nuvobea/pmat_and_mpmat/CH12F3/pmat_and_mpmat/Filtered_pmat"
# Thresholds to iterate over
THRESHOLDS="3 4 5"
# Set tolerance for comparison (0 for exact byte-by-byte match)
TOLERANCE=0

# --- Function to get file size in bytes ---
# This is necessary for exact arithmetic comparison.
get_file_size() {
    # stat -c "%s" gets the size in bytes
    stat -c "%s" "$1" 2>/dev/null || { echo "[ERROR] Could not get size for file: $1" >&2; return 1; }
}

# --- Function to get human-readable size for display ---
# This is for user output only.
get_human_size() {
    # ls -sh reports the size in a human-readable format (e.g., 71M)
    ls -sh "$1" | awk '{print $1}'
}

# --- Main Logic ---

echo "[INFO] Starting file size check iteration..."

# Iterate over sample directories
for SAMPLE_DIR in "${BASE_DIR}"/*/; do
    
    # Check if this is a directory and not an empty loop iteration
    if [ ! -d "$SAMPLE_DIR" ]; then
        continue
    fi
    
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "------------------------------------------------------------------"
    echo "[INFO] Processing sample: $SAMPLE"
    
    # Iterate over thresholds
    for THRESHOLD in $THRESHOLDS; do
        
        # Construct full file paths
        # Note: We append the directory component to BASE_DIR and then the subdirectory name to SAMPLE_DIR
        CT_BED="${SAMPLE_DIR}filtered_CT/threshold_${THRESHOLD}/filtered2_CT_${SAMPLE}_${THRESHOLD}.bed"
        GA_BED="${SAMPLE_DIR}filtered_GA/threshold_${THRESHOLD}/filtered2_GA_${SAMPLE}_${THRESHOLD}.bed"
        CTGA_BED="${SAMPLE_DIR}merged_CTGA/threshold_${THRESHOLD}/${SAMPLE}_merged_CTGA_${THRESHOLD}.bed"

        # Check for file existence and immediately report/skip if any are missing
        # We use a combined check to avoid repeating the error message for each file
        if [ ! -f "$CT_BED" ] || [ ! -f "$GA_BED" ] || [ ! -f "$CTGA_BED" ]; then
            echo "[WARNING] Skipping check for $SAMPLE at threshold $THRESHOLD: One or more required files not found."
            continue
        fi
        
        # Get file sizes in bytes (for arithmetic comparison)
        SIZE_CT=$(get_file_size "$CT_BED")
        SIZE_GA=$(get_file_size "$GA_BED")
        SIZE_CTGA=$(get_file_size "$CTGA_BED")
        
        # Calculate expected size
        EXPECTED_CTGA=$((SIZE_CT + SIZE_GA))
        
        # Calculate difference (absolute value)
        DIFFERENCE=$((SIZE_CTGA - EXPECTED_CTGA))
        if [ "$DIFFERENCE" -lt 0 ]; then
            DIFFERENCE=$((-DIFFERENCE))
        fi
        
        # Get human-readable sizes for output
        HUMAN_CT=$(get_human_size "$CT_BED")
        HUMAN_GA=$(get_human_size "$GA_BED")
        HUMAN_CTGA=$(get_human_size "$CTGA_BED")
        
        # Check if the difference is within the tolerance (which is 0 for exact match)
        if [ "$DIFFERENCE" -le "$TOLERANCE" ]; then
            echo "[INFO] Merging check PASSED for $SAMPLE (T=$THRESHOLD)."
            echo "[INFO] Size: CT ($HUMAN_CT) + GA ($HUMAN_GA) = Merged ($HUMAN_CTGA)."
        else
            echo "[ERROR] Merging check FAILED for $SAMPLE (T=$THRESHOLD)."
            echo "[ERROR] Expected total bytes: $EXPECTED_CTGA | Found bytes: $SIZE_CTGA."
            echo "[ERROR] Difference: $DIFFERENCE bytes. Check concatenation integrity."
        fi

    done
    echo "---"
done

echo "[INFO] Script finished successfully."
exit 0