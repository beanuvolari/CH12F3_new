#!/bin/bash

# nohup ./01_0SharedHotspots.sh CH12F3 3 5 0 <samples_in_order> > nohup_SharedHotspots35_merge0.out 2>&1 &

set -o errexit
set -o pipefail

##########################
#       FUNCTIONS        #
##########################

source /scratch/scripts/General_Functions.sh

# --- PATHS CONFIGURATION ---
BASE_DIR="/scratch/pmat/Filtered_pmat"
OUTPUT_BASE="/scratch/DepIndep_dataset"

# Check input parameters
if [ "$#" -lt 4 ]; then
    echo "[ERROR] Usage: $0 <cell_line> <threshold_num> <min_count> <enlargement>"
    exit 1
fi

# Input parameters
cell_line="$1"
threshold_num="$2"
min_count="$3"
enlargement="$4"
shift 4
# Remaining arguments are sample names in order
ORDERED_SAMPLES=("$@")
echo "[DEBUG] Ordered samples: ${ORDERED_SAMPLES[*]}"


if [ ! -d "$BASE_DIR" ]; then
    echo "[ERROR] La cartella radice $BASE_DIR non esiste!"
    exit 1
fi

# Global array to collect all rider file paths
all_rider_paths=()

##########################
#          MAIN          #
##########################

start_timer "SharedHotspots_merge"
echo "[INFO] Scanning for samples in: $BASE_DIR"
echo "[INFO]  Cell Line: $cell_line | Threshold: $threshold_num | MinCount: $min_count | Enlargement: $enlargement"


# 1. Use ORDERED_SAMPLES instead of scanning directories
sample_dirs=()

for sample in "${ORDERED_SAMPLES[@]}"; do
    sample_path="${BASE_DIR}/${sample}"
    if [ -d "$sample_path" ]; then
        sample_dirs+=("$sample_path")
    else
        echo "[WARNING] Directory not found for sample: $sample"
    fi
done

if [ "${#sample_dirs[@]}" -eq 0 ]; then
    echo "[ERROR] Nessuna directory valida trovata tramite ORDERED_SAMPLES!"
    exit 1
fi


# 2. Iterate over each sample directory to find the specific rider files
for sample_path in "${sample_dirs[@]}"; do
    sample_name=$(basename "$sample_path")

    target_dir="${sample_path}/merged_CTGA/threshold_${threshold_num}/${threshold_num}_rider${min_count}"
    
    # File pattern to match
    file_pattern="${target_dir}/${sample_name}_merged_CTGA_${threshold_num}_sorted-RIDER.clean*.bed"
    
    # search for matching files
    files=( $file_pattern )
    
    if [ -f "${files[0]}" ]; then
        echo "[DEBUG] Found: ${sample_name}"
        all_rider_paths+=("${files[0]}")
    else
        echo "[WARNING] File not found for sample: ${sample_name} inside $target_dir"
        true
    fi
done

# Check if any files were found
if [ "${#all_rider_paths[@]}" -eq 0 ]; then
    echo "[ERROR] No rider files found matching the criteria. Exiting."
    exit 1
fi

echo -e "[INFO] Total files collected: ${#all_rider_paths[@]}"

# Prepare output folder and file
output_dir="${OUTPUT_BASE}/Threshold_${threshold_num}/Enlargement_${enlargement}"
mkdir -p "$output_dir"
output_file="${output_dir}/${enlargement}_${cell_line}_SharedHotspots_${threshold_num}_rider${min_count}.bed"

echo "[INFO] Launching R script..."

# Call the R script
Rscript /scratch/scripts/01_1SharedHotspots_merge.R "${all_rider_paths[@]}" "$output_file" "$enlargement"

echo "[INFO] Result saved in $output_file"
end_timer "SharedHotspots_merge"
echo "[INFO] Processing completed!"

