#!/bin/bash
# nohup ./04_Merging_CTGA.sh > nohup_04Merging_CTGA.out 2>&1 &

###########################################
#              FUNCTIONS                  #
###########################################

# Source general utility functions
source /scratch/scripts/General_Functions.sh

# --- File Size Utility Functions ---
TOLERANCE=0

get_file_size() {
    stat -c "%s" "$1" 2>/dev/null || { echo "[ERROR] Could not get size for file: $1" >&2; return 1; }
}

get_human_size() {
    ls -sh "$1" | awk '{print $1}'
}

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

# Function to check if the merged file size equals the sum of its components
check_merging_integrity() {
    local CT_BED="$1"
    local GA_BED="$2"
    local MERGED_BED="$3"
    local threshold="$4"
    local filebase="$5"
    
    echo -e "\n[INFO] Starting integrity check for merged BED file (T=$threshold)."

    # Check existence again (redundant but safe)
    if [ ! -f "$CT_BED" ] || [ ! -f "$GA_BED" ] || [ ! -f "$MERGED_BED" ]; then
        echo "[WARNING] Cannot perform integrity check: one or more files are missing for T=$threshold."
        return 1
    fi

    # Get file sizes in bytes
    local SIZE_CT=$(get_file_size "$CT_BED")
    local SIZE_GA=$(get_file_size "$GA_BED")
    local SIZE_MERGED=$(get_file_size "$MERGED_BED")
    local EXPECTED_MERGED=$((SIZE_CT + SIZE_GA))
    local DIFFERENCE=$((SIZE_MERGED - EXPECTED_MERGED))
    [ "$DIFFERENCE" -lt 0 ] && DIFFERENCE=$((-DIFFERENCE))
    
    local HUMAN_CT=$(get_human_size "$CT_BED")
    local HUMAN_GA=$(get_human_size "$GA_BED")
    local HUMAN_MERGED=$(get_human_size "$MERGED_BED")
    
    # Check if the difference is within the tolerance
    if [ "$DIFFERENCE" -le "$TOLERANCE" ]; then
        echo "[INFO] Integrity check PASSED for $filebase (T=$threshold)."
        echo "[INFO] Size: CT ($HUMAN_CT) + GA ($HUMAN_GA) = Merged ($HUMAN_MERGED)."
        return 0 # Success
    else
         echo "[ERROR] Integrity check FAILED for $filebase (T=$threshold). Diff: $DIFFERENCE bytes."
         return 1 # Failure
    fi
}

# Function that merges CT and GA BED files, sorts them, and saves the output
merge_and_sort_BED_files() {
  local sample_dir="$1"
  local filebase="$2"
  local threshold="$3"
  local output_dir="$4"

  local ct_file="${sample_dir}/filtered_CT/threshold_${threshold}/filtered2_CT_${filebase}_${threshold}.bed"
  local ga_file="${sample_dir}/filtered_GA/threshold_${threshold}/filtered2_GA_${filebase}_${threshold}.bed"
  local merged_file="${output_dir}/${filebase}_merged_CTGA_${threshold}.bed"
  local merged_sorted_file="${output_dir}/${filebase}_merged_CTGA_${threshold}_sorted.bed"

  # Remove the merged file if it already exists
  rm -f "$merged_file" "$merged_sorted_file"

  echo -e "\n"
  start_timer "Merge_and_Sort_bed_file"

   if [[ -f "$ct_file" && -f "$ga_file" ]]; then
        echo "[INFO] Merging $ct_file and $ga_file"
        cat "$ct_file" "$ga_file" > "$merged_file"
        
        # Integrity check
        check_merging_integrity "$ct_file" "$ga_file" "$merged_file" "$threshold" "$filebase"

        echo "[INFO] Sorting merged BED file..."
        if sort -k1,1V -k2,2n -k3,3n "$merged_file" > "$merged_sorted_file"; then
            echo "[INFO] Sorted successfully: $merged_sorted_file"
            rm "$merged_file"
            indexing_bed "$merged_sorted_file"
        else
            echo "[ERROR] Sorting failed for $merged_file"
        fi
    else
        echo "[ERROR] Source BED files missing for threshold $threshold"
    fi

  # Check if output was generated correctly
  check_output_generation "$merged_sorted_file" 
  end_timer "Merge_and_Sort_bed_file"

}

# Function that merges variableStep WIG files from CT and GA, sorts them, and creates a BigWig
merge_variableStep_wigs() {
  local sample_dir="$1"
  local filebase="$2"
  local threshold="$3"
  local output_dir="$4"
  local chrom_sizes="$5"

  local combined_wig="${output_dir}/${filebase}_merged_CTGA_${threshold}.wig"
  local combined_bw="${output_dir}/${filebase}_merged_CTGA_${threshold}.bw"
  local temp_data=$(mktemp)
    
  rm -f "$combined_wig" "$combined_bw"
  start_timer "Merge_wig_file_and_generate_bw"

  # Loop over CT and GA wig files
  for type in CT GA; do
    local wig_file="${sample_dir}/filtered_${type}/threshold_${threshold}/filtered2_${type}_${filebase}_${threshold}.wig"
    if [ -f "$wig_file" ]; then
      echo "[INFO] Extracting data from $wig_file"

      # Extract data lines (skip variableStep headers), print: chrom pos value
      awk '
        /^variableStep chrom=/ {
          split($2,a,"="); chrom=a[2];
          next;
        }
        /^[0-9]/ {print chrom, $1, $2}
      ' "$wig_file" >> "$temp_data"
    else
      echo "[ERROR] WIG file $wig_file not found!"
    fi
  done

  # Sort data by chromosomal order and position
     if [ -s "$temp_data" ]; then
        echo "[INFO] Sorting and formatting WIG data..."
        local temp_sorted=$(mktemp)
        sort -k1,1V -k2,2n "$temp_data" > "$temp_sorted"

        awk 'BEGIN {last_chrom=""} {
            if ($1 != last_chrom) {
                print "variableStep chrom=" $1 " span=1";
                last_chrom = $1;
            }
            print $2 "\t" $3;
        }' "$temp_sorted" > "$combined_wig"

        check_output_generation "$combined_wig"  # Check if output was generated correctly
   
        # Create BigWig file if chromosome sizes file exists
        if [[ -f "$chrom_sizes" && -f "$combined_wig" ]]; then
            if wigToBigWig "$combined_wig" "$chrom_sizes" "$combined_bw"; then
                echo "[INFO] BigWig file created: $combined_bw"
            else
                echo "[ERROR] wigToBigWig conversion failed for $combined_wig."
            fi
        else
          echo "[ERROR] Chromosome sizes file not found at $chrom_sizes — skipping BigWig conversion"
        fi
        rm "$temp_data" "$temp_sorted"
    else
        echo "[ERROR] No data found to merge for WIG files."
    fi

    end_timer "Merge_wig_file_and_generate_bw"
}

###########################################
#                 MAIN                    #
###########################################

# Define directories and files
pmat_dir="scratch/pmat/Filtered_pmat"
chrom_sizes="/scratch/genome/mm9.chrom.sizes"

# Check if directory exists
check_directory_exists "$pmat_dir"

# Loop over each sample directory inside pmat_dir
for sample_dir in "${pmat_dir}"/*; do
  if [ -d "$sample_dir" ]; then
    echo -e "[INFO] Entered directory: $sample_dir\n"

    filebase=$(basename "$sample_dir")
    echo -e "\n--------------------------------------------------------------------------------------\n"
    echo -e "\n[INFO] Starting sample: $filebase"

    dir_merged="${sample_dir}/merged_CTGA"
    check_directory_exists "$dir_merged"

  # Threshold detection
    if [ "$#" -gt 0 ]; then
      # Arguments provided: parse them as thresholds
      input_thresholds="$*"
      read -r -a thresholds <<< "$(parse_thresholds "$input_thresholds")"
      echo "[INFO] Thresholds specified as arguments: ${thresholds[*]}"
    else
      # No arguments: detect threshold directories automatically
      # Find all threshold directories inside filtered_CT, sorted numerically
      echo "[INFO] No threshold arguments provided — scanning directories"
      mapfile -t threshold_dirs < <(find "${sample_dir}/filtered_CT" -maxdepth 1 -type d -name "threshold_*" -exec basename {} \; | sort -V)
      thresholds=()
      for dir in "${threshold_dirs[@]}"; do
        thresholds+=("${dir#threshold_}")
      done
    fi
    echo -e "[INFO] Thresholds to process: ${thresholds[*]}"

    # Iterate over thresholds for processing
    for threshold in "${thresholds[@]}"; do
      echo -e "\n-------------------------"
      echo -e "\n[INFO] Processing threshold: $threshold"
      
      # Create directory for merged files at this threshold
      dir_merged_threshold="${dir_merged}/threshold_${threshold}"
      check_directory_exists "$dir_merged_threshold"

      merge_and_sort_BED_files "$sample_dir" "$filebase" "$threshold" "$dir_merged_threshold"
      merge_variableStep_wigs "$sample_dir" "$filebase" "$threshold" "$dir_merged_threshold" "$chrom_sizes"
    done

    echo -e "[INFO] Analysis completed for sample: $filebase"
    echo -e "\n--------------------------------------------------------------------------------------\n"
  fi
done

echo "[INFO] All tasks finished successfully."
