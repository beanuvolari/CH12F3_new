#!/bin/bash
# nohup ./05_Index_forIGVbrowser.sh > nohup_05_IndexIGVbrowser.out 2>&1 &

set -euo pipefail



###########################################
#              FUNCTIONS                  #
###########################################
# source /scratch/scripts/General_Functions.sh 
source /scratch/pmat/Scripts/General_Functions.sh 

indexing_bed() {
  local input_bed="$1"
  local gz_file="${input_bed}.gz"
  local index_file="${gz_file}.tbi"

  if [[ -f "$input_bed" ]]; then
    echo "[INFO] Processing BED file: $input_bed"

    # if file '.gz' already exist, skip it
    if [[ -f "$gz_file" && -f "$index_file" ]]; then
      echo "[INFO] Already compressed and indexed: $gz_file"
      return 0
    fi

    echo "[INFO] Compressing with bgzip..."
    bgzip -c -@ "$(nproc)" "$input_bed" > "$gz_file"

    echo "[INFO] Creating tabix index..."
    tabix -p bed "$gz_file"

    if [[ -f "$index_file" ]]; then
      echo "[INFO] Index created successfully: $index_file"
    else
      echo "[ERROR] Failed to create index for: $gz_file"
    fi
  else
    echo "[ERROR] File not found: $input_bed"
  fi
}



###########################################
#                 MAIN                    #
###########################################

# base_dir="/scratch/pmat/Filtered_pmat"
base_dir="/scratch/pmat/pmat_and_mpmat/Filtered_pmat"

start_timer "Sorting_bed"

echo -e "\n[INFO] Starting RIDER batch processing in: $base_dir"

for sample_dir in "$base_dir"/*; do
    merged_dir="$sample_dir/merged_CTGA"

    # salta se merged_CTGA non esiste
    [[ -d "$merged_dir" ]] || continue

    for threshold_dir in "$merged_dir"/threshold_*; do
        [[ -d "$threshold_dir" ]] || continue

        for bed_file in "$threshold_dir"/*sorted.bed; do
            [[ -f "$bed_file" ]] || continue

            echo -e "\n--------------------------------------"
            echo "[INFO] Found BED file: $bed_file"
            indexing_bed "$bed_file"
        done
    done
done


end_timer "Sorting_bed"

echo -e "\n[INFO] Process completed successfully"