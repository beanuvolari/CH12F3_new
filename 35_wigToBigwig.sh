#!/bin/bash
# nohup ./35_wigToBigwig.sh > nohup_35_wigToBigwig.out 2>&1 &
###########################################
#              FUNCTIONS                  #
###########################################

# Source general utility functions
source /scratch/pmat/General_Functions.sh  

# Function that sorts BED files and saves the output
sort_BED_files() {
  local filebase="$1"
  local threshold="$2"
  local working_dir="$3"
  local chrom_sizes="$4"
  local type="$5"

  # Output file
  local sorted_file="${working_dir}/filtered2_${type}_${filebase}_${threshold}_sorted.bed"

  # Remove the sorted file if it already exists
  [ -f "$sorted_file" ] && rm "$sorted_file"

  echo -e "\n"
  start_timer "Sort_bed_file"

  # Input file
  local input_file="${working_dir}/filtered2_${type}_${filebase}_${threshold}.bed"
    if [ -f "$input_file" ]; then
    echo "[INFO] Processing $input_file"
    
    # Sort the BED file by chromosome and start position
    if bedtools sort -faidx "$chrom_sizes" -i "$input_file" > "$sorted_file"; then
      echo "[INFO] Sorted successfully: $sorted_file"
    else
      echo "[ERROR] Sorting failed for $input_file"
    fi
  else
    echo -e "[ERROR] File $input_file not found!"
  fi
  
  check_output_generation "$sorted_file"  # Check if output was generated correctly
  end_timer "Sort_bed_file"
}

# Function creates a BigWig
generate_bw() {
  local filebase="$1"
  local threshold="$2"
  local working_dir="$3"
  local chrom_sizes="$4"
  local type="$5"

  # Output file
  local bw_file="${working_dir}/filtered2_${type}_${filebase}_${threshold}.bw"

  # Remove the BW file if it already exists
  [ -f "$bw_file" ] && rm "$bw_file"
  
  start_timer "Generate_bw"

  # Input file
    local wig_file="${working_dir}/filtered2_${type}_${filebase}_${threshold}.wig"
    if [ -f "$wig_file" ]; then
      echo "[INFO] Processing $wig_file"

      # Create BigWig file if chromosome sizes file exists
      if [ -f "$chrom_sizes" ]; then
        if wigToBigWig "$wig_file" "$chrom_sizes" "$bw_file"; then
          echo "[INFO] BigWig file created: $bw_file"
        else
          echo "[ERROR] wigToBigWig conversion failed for $bw_file"
        fi
      else
        echo "[ERROR] Chromosome sizes file not found at $chrom_sizes — skipping BigWig conversion"
      fi
    else
      echo -e "[ERROR] File $wig_file not found!"
    fi
  end_timer "Generate_bw"
}
###########################################
#                 MAIN                    #
###########################################

# Directory containing pmat data
pmat_dir="/scratch/pmat/Prova_Pipeline"

# Check if directory exists
check_directory_exists "$pmat_dir"

# Define chrom_size file
chrom_sizes="/scratch/pmat/Controls/mm9.chrom.sizes"

# Folders to be skipped
skip_samples=("Dep_Indep_dataset" "Fabio_images" "prove_AIDDependency" "Threshold_3" "Threshold_4" "Threshold_5")


# Loop over each sample directory inside pmat_dir
for sample_dir in "${pmat_dir}"/*; do
  if [ -d "$sample_dir" ]; then
  filebase=$(basename "$sample_dir")
    
     # If the sample is in the skip list, continue to the next iteration
    if [[ "${skip_samples[@]}" =~ "${filebase}" ]]; then
      echo "[INFO] Skipping sample: $filebase"
      continue
    fi

    echo -e "[INFO] Entered directory: $sample_dir\n"

    echo -e "\n--------------------------------------------------------------------------------------\n"
    echo -e "\n[INFO] Analysis started for sample: $filebase"

    for type in CT GA; do
        dir_filtered="${sample_dir}/filtered_${type}"
        check_directory_exists "$dir_filtered"  # Create output directory if missing

      # Find all threshold directories inside filtered_CT, sorted numerically
      mapfile -t threshold_dirs < <(find "${sample_dir}/filtered_${type}" -type d -name "threshold_*" \
      -exec basename {} \; | sort -V)

      for threshold_dir in "${threshold_dirs[@]}"; do
        threshold="${threshold_dir#threshold_}"  # Extract threshold number
        echo -e "\n-------------------------"
        echo -e "\n[INFO] Processing threshold: $threshold"
        
        # Create directory for merged files at this threshold
        dir_filtered_threshold="${dir_filtered}/threshold_${threshold}"
        check_directory_exists "$dir_filtered_threshold"

        sort_BED_files "$filebase" "$threshold" "$dir_filtered_threshold" "$chrom_sizes" "$type"
        generate_bw "$filebase" "$threshold" "$dir_filtered_threshold" "$chrom_sizes" "$type"
      done

      echo -e "[INFO] Analysis completed for sample: $filebase"
      echo -e "\n--------------------------------------------------------------------------------------\n"
    done
  fi
done

echo "[INFO] All analyses have been completed!"
