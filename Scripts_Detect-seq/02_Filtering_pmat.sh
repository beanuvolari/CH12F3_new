#!/bin/bash

###########################################
# This script requires as input: PMAT files.
#  Part 1:
#   - filters PMAT files to extract CT and GA mutations;
#   - generates WIG files for CT and GA mutations.
#  Part 2:
#   - applies thresholds to the CT and GA mutations;
#   - generates filtered PMAT, WIG, and BED files for each threshold.
#  It can be run interactively or with command line arguments for thresholds.
#  The script is designed to be run in a batch mode, and it can handle multiple samples.
###########################################

# Usage: nohup ./01_Filtering_pmat.sh <threshold1> <threshold2> ... > nohup_01_Filtering.out 2>&1 &
# Example: nohup ./01_Filtering_pmat.sh 3 4 5 > nohup_01_Filtering2.out 2>&1 &
#          nohup ./01_Filtering_pmat.sh --part1 > nohup_01_Filtering_part1.out 2>&1 &       -> Only filtering for CT and GA
#          nohup ./01_Filtering_pmat.sh --part2 4 5 > nohup_01_Filtering_part2.out 2>&1 &   -> Only filtering by thresholds

###########################################
#           GLOBAL VARIABLES              #
###########################################

declare -gA chrom_written  # Global associative array for tracking chromosomes in WIG files

###########################################
#              FUNCTIONS                  #
###########################################

source /scratch/scripts/General_Functions.sh

# Function to write WIG files from an associative array
append_to_wig() {
    local wig_file="$1"     # WIG file output
    local chrom="$2"        # Chromosome name (e.g.: chr1)
    local pos="$3"          # Genomic position (e.g.: 123456)
    local value="$4"        # Coverage

     # If the wig file does not already have the header with the cromosome name, add it
    if [[ -z "${chrom_written[$wig_file,$chrom]}" ]]; then
        echo "variableStep chrom=$chrom span=1" >> "$wig_file"
        chrom_written[$wig_file,$chrom]=1
    fi

    # Add the row with position and value
    echo "$pos $value" >> "$wig_file"
}


# Function: Extract CT and GA mutations from PMAT file
filter_CT_GA_mutations() {
    local input_file="$1" # Input file need a PMAT format
    local output_dir="$2"
    local sample="$3" # Sample name extracted from the input file
    
    # Reset the global array for tracking chromosomes in WIG files
    chrom_written=() 
    
    # Check if the input file exists
    check_file_exists "$input_file"
    
    # Check if the output directory exists, if not create it
    check_directory_exists "$output_dir"

    # Create output directories for the sample
    local sample_dir="$output_dir/$sample"
    mkdir -p "$sample_dir/filtered_CT" "$sample_dir/filtered_GA"

    # Clean files if they already exist
    for f in \
    "${sample_dir}/filtered_CT/filtered_CT_${sample}.pmat" \
    "${sample_dir}/filtered_CT/filtered_CT_${sample}.wig" \
    "${sample_dir}/filtered_GA/filtered_GA_${sample}.pmat" \
    "${sample_dir}/filtered_GA/filtered_GA_${sample}.wig" 
    do
         > "$f"
    done

    # Starting the filtering process
    echo -e "\n[INFO] Starting CT/GA extraction from: $input_file"
    IFS=$'\t' grep -E 'CT|GA' "$input_file" | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15
    do
        local output_pmat="$sample_dir/filtered_${col9}/filtered_${col9}_${sample}.pmat"
        local output_wig="$sample_dir/filtered_${col9}/filtered_${col9}_${sample}.wig"

        # Append line to .pmat file
        printf "%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$col1" "$col2" "$col5" "$col6" "$col7" "$col8" "$col9" "$col12" "$col13" "$col14" >> "$output_pmat"

         # Determinate the CT GA value in column 9 for WIG file
        local wig_value
        if [[ $col9 == "CT" ]]; then
            wig_value="$col13"
        elif [[ $col9 == "GA" ]]; then
            wig_value="-$col13"
        fi

        # Use the append_to wig function to correctly write into WIG file
        append_to_wig "$output_wig" "$col1" "$col2" "$wig_value" 
    done
    
    # Check if the output files have been created
    check_output_generation "$sample_dir/filtered_CT/filtered_CT_${sample}.pmat"
    check_output_generation "$sample_dir/filtered_CT/filtered_CT_${sample}.wig"
    check_output_generation "$sample_dir/filtered_GA/filtered_GA_${sample}.pmat"
    check_output_generation "$sample_dir/filtered_GA/filtered_GA_${sample}.wig"

    echo -e "[INFO] Filtering CT_GA mutations completed for: $sample"
}

# Function: Apply threshold filtering to mutations
filter_by_threshold() {
    local input_file="$1"
    local output_dir="$2"
    local threshold="$3"
    local sample="$4"
    local mut_type="$5"  # Mutation type (CT or GA)
    
    # Reset the global array for tracking chromosomes in WIG files
    chrom_written=()  

    # Check if the input file exists
    check_file_exists "$input_file"

    # Clean files if they already exist
    for f in \
        "${output_dir}/filtered_${mut_type}/threshold_${threshold}/filtered2_${mut_type}_${sample}_${threshold}.pmat" \
        "${output_dir}/filtered_${mut_type}/threshold_${threshold}/filtered2_${mut_type}_${sample}_${threshold}.wig" \
        "${output_dir}/filtered_${mut_type}/threshold_${threshold}/filtered2_${mut_type}_${sample}_${threshold}.bed"
    do
        > "$f"
    done

    # Define the output directory for the threshold
    local threshold_dir="$output_dir/filtered_${mut_type}/threshold_${threshold}"
        check_directory_exists "$threshold_dir"

    # Start the filtering process
    echo -e "[INFO] Starting threshold filtering for: $sample $mut_type with threshold = $threshold"
    while IFS=$'\t' read -r firstField col3 col4 col5 col6 col7 col8 col9 col10
    do
        IFS='_' read -r col1 col2 <<< "$firstField"

        local output2_pmat="$threshold_dir/filtered2_${col7}_${sample}_${threshold}.pmat"
        local output2_wig="$threshold_dir/filtered2_${col7}_${sample}_${threshold}.wig"
        local output2_bed="$threshold_dir/filtered2_${col7}_${sample}_${threshold}.bed"

        # Apply threshold filtering
        if [[ $col9 -ge $threshold ]]; then
            printf "%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$col1" "$col2" "$col3" "$col4" "$col5" "$col6" "$col7" "$col8" "$col9" "$col10" >> "$output2_pmat"

            local wig_value
            if [[ $col7 == "CT" ]]; then
                wig_value="$col9"
            elif [[ $col7 == "GA" ]]; then
                wig_value="-$col9"
            fi

            append_to_wig "$output2_wig" "$col1" "$col2" "$wig_value"

            # Generate BED entries for each count
            for ((i = 0; i < col9; i++)); do
                uniqueName=$i
                randomNumber=$((RANDOM % 31 + 30))
                strand=$( (( RANDOM % 2 == 0 )) && echo "+" || echo "-")
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$col1" "$col2" "$col2" "$uniqueName" "$randomNumber" "$strand" >> "$output2_bed"
            done
        fi

    done < "$input_file"
    
    # Check if the output files have been created
    local threshold_dir="$output_dir/filtered_${mut_type}/threshold_${threshold}"
    local output2_pmat="$threshold_dir/filtered2_${mut_type}_${sample}_${threshold}.pmat"
    local output2_wig="$threshold_dir/filtered2_${mut_type}_${sample}_${threshold}.wig"
    local output2_bed="$threshold_dir/filtered2_${mut_type}_${sample}_${threshold}.bed"

    check_output_generation "$output2_pmat"
    check_output_generation "$output2_wig"
    check_output_generation "$output2_bed"

    echo -e "[INFO] Filtering by threshold = $threshold completed for: $sample $mut_type"
}

# Function: Process threshold filtering for all samples
process_thresholds_for_samples() {
    local output_dir="$1"
    local -n thresholds_ref=$2  # reference to array of thresholds

    # Check if the output directory exists, if not create it
    check_directory_exists "$output_dir"

    # Iterate over all sample directories in the base directory
    for sample_dir in "$output_dir"/*; do
        if [ -d "$sample_dir" ]; then
            local sample_name=$(basename "$sample_dir")
            echo -e "\n-------------------------"
            echo -e "\n[INFO] Starting analysis for: $sample_name"

            # Record start time
            start_timer "$sample_name"

            # Iterate over each threshold
            for threshold in "${thresholds_ref[@]}"; do
             echo -e "\n-----"
                echo -e "\n[INFO] Applying threshold: $threshold"

                # Generate the outputs with CT and GA mutations
                for mut_type in CT GA; do
                    echo -e "\nProcessing type: $mut_type"
                    start_timer "$mut_type"
                    file_path="$sample_dir/filtered_${mut_type}/filtered_${mut_type}_${sample_name}.pmat"
                    if [ -f "$file_path" ]; then
                        filter_by_threshold "$file_path" "$sample_dir" "$threshold" "$sample_name" "$mut_type"
                    else
                        echo "[WARN] File not found: $file_path"
                    fi
                    end_timer "$mut_type"
                    echo -e "\n"                  
                done
            done
            
            # Record end time
            echo -e "\n"
            end_timer "$sample_name"
        fi
    done
}

###########################################
#                 MAIN                    #
###########################################

# Default run mode
run_part1=true
run_part2=true

# Check for mode flags
case "$1" in
    --part1)
        run_part2=false
        shift
        ;;
    --part2)
        run_part1=false
        shift
        ;;
    --all)
        shift
        ;;
esac


# Define input/output directories
pmat_dir="/scratch/pmat"
output_dir="${pmat_dir}/Filtered_pmat"

# Record the start time
echo -e "\n--------------------------------------------------------------------------------------\n"
start_timer "Analysis_Filtering_pmat"

# Check if the input directory exists
check_directory_exists "$output_dir"

### Part 1: Extract CT/GA mutations from all .pmat files
if [ "$run_part1" = true ]; then
    echo -e "\n-------------------------"
    echo -e "\n\nPart 1: Extracting CT and GA mutations from PMAT files"
    echo -e "\n-----------------------------------------------------------------\n"

    # Iterate over all .pmat files in the input directory
    for input_pmat in "${pmat_dir}"/*.pmat; do
        sample=$(basename "$input_pmat" .pmat)
            # Remove _merge_CLEAN or _CLEAN suffixes from sample name
            sample=${sample/_merge_CLEAN/}
            sample=${sample/_CLEAN/}
            echo -e "\n[INFO] Sample name extracted: $sample"

        # Record the start time
        start_timer "$sample"

        # Start filtering CT and GA mutations
        filter_CT_GA_mutations "$input_pmat" "$output_dir" "$sample"

        # Record the end time
        echo -e "\n"
        end_timer "$sample"
        echo -e "-------------------------\n"
    done
fi


### Part 2: Handle interactive or argument-based threshold input
if [ "$run_part2" = true ]; then   
    echo -e "\n\nPart 2: Applying thresholds to CT and GA mutations"
    echo -e "\n-----------------------------------------------------------------\n"


    # Check if thresholds are provided as arguments
    # If no arguments are provided, prompt for thresholds
    if [ "$#" -lt 1 ]; then
        read -p "Enter threshold list or range (e.g. '3 5 7' or '3-6'): " input_thresholds
    else
        input_thresholds="$*"
    fi

    # Parse the thresholds
    thresholds=($(parse_thresholds "$input_thresholds"))

    echo -e "\n[INFO] Thresholds to be applied: ${thresholds[@]}"
    echo -e "\n[INFO] Starting threshold-based filtering for all samples...\n"

    # Start threshold-based filtering for all samples
    process_thresholds_for_samples "$output_dir" thresholds
    echo -e "\n-----------------------------------------------------------------\n"
fi

# Record the end time
end_timer "Analysis_Filtering_pmat"
echo -e "\n--------------------------------------------------------------------------------------\n"

echo -e "End of the process.\n"