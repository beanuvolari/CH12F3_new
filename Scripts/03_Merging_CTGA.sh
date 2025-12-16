#!/bin/bash
# nohup ./03_Merging_CTGA.sh > nohup_03_Merging_CTGA3.out 2>&1 &
###########################################
#              FUNCTIONS                  #
###########################################

# Source general utility functions
source /scratch/scripts/General_Functions.sh  

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


# --- File Size Utility Functions ---
# [INFO] Setting tolerance for byte size check. 0 for exact match.
TOLERANCE=0

get_file_size() {
    # stat -c "%s" gets the size in bytes
    stat -c "%s" "$1" 2>/dev/null || { echo "[ERROR] Could not get size for file: $1" >&2; return 1; }
}

get_human_size() {
    # ls -sh reports the size in a human-readable format (e.g., 71M)
    ls -sh "$1" | awk '{print $1}'
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
    SIZE_CT=$(get_file_size "$CT_BED")
    SIZE_GA=$(get_file_size "$GA_BED")
    SIZE_MERGED=$(get_file_size "$MERGED_BED")
    
    # Calculate expected size
    EXPECTED_MERGED=$((SIZE_CT + SIZE_GA))
    
    # Calculate difference (absolute value)
    DIFFERENCE=$((SIZE_MERGED - EXPECTED_MERGED))
    if [ "$DIFFERENCE" -lt 0 ]; then
        DIFFERENCE=$((-DIFFERENCE))
    fi
    
    # Get human-readable sizes for output
    HUMAN_CT=$(get_human_size "$CT_BED")
    HUMAN_GA=$(get_human_size "$GA_BED")
    HUMAN_MERGED=$(get_human_size "$MERGED_BED")
    
    # Check if the difference is within the tolerance
    if [ "$DIFFERENCE" -le "$TOLERANCE" ]; then
        echo "[INFO] Integrity check PASSED for $filebase (T=$threshold)."
        echo "[INFO] Size: CT ($HUMAN_CT) + GA ($HUMAN_GA) = Merged ($HUMAN_MERGED)."
        return 0 # Success
    else
        echo "[ERROR] Integrity check FAILED for $filebase (T=$threshold)."
        echo "[ERROR] Expected total bytes: $EXPECTED_MERGED | Found bytes: $SIZE_MERGED."
        echo "[ERROR] Difference: $DIFFERENCE bytes. Concatenation integrity failed."
        return 1 # Failure
    fi
}

# Function that merges CT and GA BED files, sorts them, and saves the output
merge_and_sort_BED_files() {
  local sample_dir="$1"
  local filebase="$2"
  local threshold="$3"
  local output_dir="$4"
  local chrom_sizes="$5"

  local merged_file="${output_dir}/${filebase}_merged_CTGA_${threshold}.bed"
  local merged_sorted_file="${output_dir}/${filebase}_merged_CTGA_${threshold}_sorted.bed"

  # Remove the merged file if it already exists
  [ -f "$merged_file" ] && rm "$merged_file"
  [ -f "$merged_sorted_file" ] && rm "$merged_sorted_file"

  echo -e "\n"
  start_timer "Merge_and_Sort_bed_file"

  # Loop over CT and GA file types
  for type in CT GA; do
    local bed_file="${sample_dir}/filtered_${type}/threshold_${threshold}/filtered2_${type}_${filebase}_${threshold}.bed"
    if [ -f "$bed_file" ]; then
      echo -e "[INFO] Adding $bed_file to $merged_file"
      cat "$bed_file" >> "$merged_file"
    else
      echo -e "[ERROR] File $bed_file not found!"
    fi
  done

  # Perform integrity check on the UNSORTED file before sorting
    check_merging_integrity "$CT_BED_IN" "$GA_BED_IN" "$merged_file_unsorted" "$threshold" "$filebase"

  # Sort the merged BED file by chromosome and start position
  if [ -f "$merged_file_unsorted" ]; then
    echo "[INFO] Sorting merged BED file: $merged_file_unsorted"
    
    if sort -k1,1V -k2,2n -k3,3n -k4,4n "$merged_file_unsorted" > "$merged_sorted_file"; then
      echo "[INFO] Sorted successfully: $merged_sorted_file"

      # Remove the unsorted merged file to save space
      rm "$merged_file_unsorted"

      # Index the final sorted file. This ensures the original sorted file remains.
      indexing_bed "$merged_sorted_file"

    else
      echo "[ERROR] Sorting failed for $merged_file_unsorted"
    fi
  else
    echo "[ERROR] Merged BED file $merged_file_unsorted not found for sorting!"
  fi

  
  check_output_generation "$merged_sorted_file"  # Check if output was generated correctly
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

  # Remove existing combined wig and bw files if they exist
  [ -f "$combined_wig" ] && rm "$combined_wig"
  [ -f "$combined_bw" ] && rm "$combined_bw"

  # Temporary file to collect chrom, position, and value data
  local temp_data
  temp_data=$(mktemp)
  
  start_timer "Merge_wig_file_and_generate_bw"

  # Loop over CT and GA wig files
  for type in CT GA; do
    local wig_file="${sample_dir}/filtered_${type}/threshold_${threshold}/filtered2_${type}_${filebase}_${threshold}.wig"
    if [ -f "$wig_file" ]; then
      echo "[INFO] Processing $wig_file"

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
  temp_sorted=$(mktemp)
  sort -k1,1V -k2,2n "$temp_data" > "$temp_sorted"

  # Write merged wig file, grouping by chromosome and adding variableStep headers
  local last_chrom=""
  while read -r chrom pos val; do
    if [[ "$chrom" != "$last_chrom" ]]; then
      echo "variableStep chrom=$chrom span=1" >> "$combined_wig"
      last_chrom="$chrom"
    fi
    echo -e "${pos}\t${val}" >> "$combined_wig"
  done < "$temp_sorted"

    #cat "${temp_data}_sorted"

  # Clean up temporary files
  rm "$temp_data" "$temp_sorted"

  check_output_generation "$combined_wig"  # Check if output was generated correctly

  # Create BigWig file if chromosome sizes file exists
  if [ -f "$chrom_sizes" ]; then
    if wigToBigWig "$combined_wig" "$chrom_sizes" "$combined_bw"; then
      echo "[INFO] BigWig file created: $combined_bw"
    else
      echo "[ERROR] wigToBigWig conversion failed for $combined_wig"
    fi
  else
    echo "[ERROR] Chromosome sizes file not found at $chrom_sizes — skipping BigWig conversion"
  fi

    end_timer "Merge_wig_file_and_generate_bw"
}


###########################################
#                 MAIN                    #
###########################################

# Directory containing pmat data
pmat_dir="scratch/pmat/Filtered_pmat"

# Check if directory exists
check_directory_exists "$pmat_dir"

# Define chrom_size file
chrom_sizes="/scratch/genome/mm9.chrom.sizes"


# Loop over each sample directory inside pmat_dir
for sample_dir in "${pmat_dir}"/*; do
  if [ -d "$sample_dir" ]; then
    echo -e "[INFO] Entered directory: $sample_dir\n"

    filebase=$(basename "$sample_dir")
    echo -e "\n--------------------------------------------------------------------------------------\n"
    echo -e "\n[INFO] Analysis started for sample: $filebase"

    dir_merged="${sample_dir}/merged_CTGA"
    check_directory_exists "$dir_merged"  # Create output directory if missing

    if [ "$#" -gt 0 ]; then
      # Arguments provided: parse them as thresholds
      input_thresholds="$*"
      read -r -a thresholds <<< "$(parse_thresholds "$input_thresholds")"
      echo "[INFO] Thresholds specified as arguments: ${thresholds[*]}"
    else
      # No arguments: detect threshold directories automatically
      # Find all threshold directories inside filtered_CT, sorted numerically
      echo "[INFO] No threshold arguments provided — scanning directories"
      mapfile -t threshold_dirs < <(find "${sample_dir}/filtered_CT" -type d -name "threshold_*" -exec basename {} \; | sort -V)
      thresholds=()
      for dir in "${threshold_dirs[@]}"; do
        thresholds+=("${dir#threshold_}")
      done
    fi


    # Iterate over thresholds for processing
    for threshold in "${thresholds[@]}"; do
      echo -e "\n-------------------------"
      echo -e "\n[INFO] Processing threshold: $threshold"
      
      # Create directory for merged files at this threshold
      dir_merged_threshold="${dir_merged}/threshold_${threshold}"
      check_directory_exists "$dir_merged_threshold"

      merge_and_sort_BED_files "$sample_dir" "$filebase" "$threshold" "$dir_merged_threshold" "$chrom_sizes"
      merge_variableStep_wigs "$sample_dir" "$filebase" "$threshold" "$dir_merged_threshold" "$chrom_sizes"
    done

    echo -e "[INFO] Analysis completed for sample: $filebase"
    echo -e "\n--------------------------------------------------------------------------------------\n"
  fi
done

echo "[INFO] All analyses have been completed!"
