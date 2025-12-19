#!/bin/bash

# nohup ./04_Rider_processing.sh --single 5 3 4 5 > nohup_04_Rider_processing_single.out 2>&1 &
# nohup ./04_Rider_processing.sh --merge 5 3-6 > nohup_04_Rider_processing_merge.out 2>&1 &

###########################################
#              FUNCTIONS                  #
###########################################

# Source general utility functions
source /scratch/scripts/General_Functions.sh  

# Function: runs RIDER on a .bed file and renames the resulting "results*" folder
run_rider_and_rename() {
  local bed_file="$1"
  local output_dir="$2"
  local threshold="$3"
  local min_count="$4"

  echo -e "\n[INFO] Running RIDER on: $bed_file"
  java -jar "$RIDER_JAR" -e_file "$bed_file" -genome "$GENOME" -filter sb=n -outDir "$output_dir" -filter min_count="$min_count"

  local results_dir
  results_dir=$(find "$output_dir" -maxdepth 1 -type d -name "results*")
  if [ -n "$results_dir" ]; then
    local new_name="${output_dir}/${threshold}_rider5"
    echo "[INFO] Renaming $results_dir to $new_name"
    mv "$results_dir" "$new_name"
  else
    echo "[WARNING] No 'results*' directory found in $output_dir"
  fi
}


# Mode: --single
run_single_mode() {
  local min_count="$1"
  shift
    local thresholds=("$@")

  start_timer "Single_mode_analysis"

  echo -e "\n[INFO] Running in SINGLE mode with min_count = $min_count"
  echo -e "[INFO] Thresholds to analyze: ${thresholds[*]}\n"

  for sample_dir in "$BASE_DIR"/*; do
    [ -d "$sample_dir" ] || continue
    sample_name=$(basename "$sample_dir")
    start_timer "$sample_name"
    echo "[INFO] Processing sample: $sample_name"

    for filter in "filtered_CT" "filtered_GA"; do
      for threshold in "${thresholds[@]}"; do
        threshold_dir="$sample_dir/$filter/threshold_${threshold}"
        [ -d "$threshold_dir" ] || {
          echo "[WARNING] Directory not found: $threshold_dir"
          continue
        }

        bed_file="${threshold_dir}/filtered2_${filter#filtered_}_${sample_name}_${threshold}.bed"
        if [ -f "$bed_file" ]; then
          run_rider_and_rename "$bed_file" "$threshold_dir" "$threshold" "$min_count"
        else
          echo "[WARNING] BED file not found: $bed_file"
        fi
      done
    done
    
    end_timer "$sample_name"

  done
  
  end_timer "Single_mode_analysis"
}

# Mode: --merge
run_merge_mode() {
  local min_count="$1"
  shift
    local thresholds=("$@")

  echo -e "\n[INFO] Running in MERGE mode with min_count = $min_count"
  echo -e "[INFO] Thresholds to analyze: ${thresholds[*]}\n"

  for sample_dir in "$BASE_DIR"/*; do
    [ -d "$sample_dir" ] || continue
    sample_name=$(basename "$sample_dir")
    start_timer "$sample_name"
    echo -e "[INFO] Processing sample: $sample_name\n"

    merged_dir="$sample_dir/merged_CTGA"
    [ -d "$merged_dir" ] || {
      echo "[WARNING] Directory not found: $merged_dir"
      continue
    }

    for threshold in "${thresholds[@]}"; do
    threshold_dir="${merged_dir}/threshold_${threshold}"
        [ -d "$threshold_dir" ] || {
          echo "[WARNING] Directory not found: $threshold_dir"
          continue
        }
      bed_file="${threshold_dir}/${sample_name}_merged_CTGA_${threshold}_sorted.bed"
      if [ -f "$bed_file" ]; then
        run_rider_and_rename "$bed_file" "$threshold_dir" "$threshold" "$min_count"
      else
        echo "[WARNING] Merged BED file not found: $bed_file"
      fi
    done

    end_timer "$sample_name"

  done
}

###########################################
#                 MAIN                    #
###########################################

RIDER_JAR="/scratch/rider/RIDER_v0.2.jar"
GENOME="mm9"
BASE_DIR="/scratch/pmat/Filtered_pmat"

# Require at least 2 arguments: mode and min_count
if [ $# -lt 2 ]; then
  echo "[ERROR] Not enough arguments."
  echo "[INFO] Usage: $0 [--single | --merge] <min_count> [thresholds | threshold_range]"
  echo "[INFO] Example: $0 --single 5 3-6"
  echo "[INFO] Example: $0 --merge 10 3 5 7"
  exit 1
fi

MODE="$1"
MIN_COUNT="$2"
shift 2

# Default to range 3-10 if not specified
#THRESHOLDS="${*:-3-10}"
# Default to range 3-10 if not specified
THRESHOLD_STRING="${*:-3-10}"
THRESHOLDS_ARRAY=($(parse_thresholds "$THRESHOLD_STRING"))


case "$MODE" in
  --single)
    run_single_mode "$MIN_COUNT" "${THRESHOLDS_ARRAY[@]}"
    ;;
  --merge)
    run_merge_mode "$MIN_COUNT" "${THRESHOLDS_ARRAY[@]}"
    ;;
  *)
    echo "[ERROR] Invalid mode: $MODE"
    echo "[INFO] Usage: $0 [--single | --merge] <min_count> [thresholds | threshold_range]"
    exit 1
    ;;
esac

echo -e "\nProcess completed"