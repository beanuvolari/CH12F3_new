#!/bin/bash

# Usage:
# nohup ./03_AID_DepIndephotspots_split.sh CH12F3 3 5 0 AID_KO_CIT > nohup_03_AID_DepIndephotspots_split350.out 2>&1 &
# parameters: <cell_line> <threshold_num> <min_count> <enlargement> <sample_name>

#########################
#       FUNCTIONS       #
#########################

source /scratch/scripts/General_Functions.sh

# Function to check if required arguments are provided
check_args() {
    if [ "$#" -lt 5 ]; then
        echo "Error: Missing arguments."
        echo "Usage: $0 <cell_line> <threshold_num> <min_count> <enlargement> <sample_name>"
        exit 1
    fi
}

# Function to construct file paths based on parameters
build_paths() {
    cell_line="$1"
    threshold_num="$2"
    min_count="$3"
    enlargement="$4"
    target_sample="$5"

    # Define base directory and input file paths
    base_dir="/scratch/DepIndep_dataset/Threshold_${threshold_num}/Enlargement_${enlargement}"
    input_file="${base_dir}/${enlargement}_${cell_line}_SharedHotspots_${threshold_num}_rider${min_count}.bed"
    check_file_exists "$input_file"

    # Define output directory and output file paths
    outdir_dep="${base_dir}/Dependent_hotspots"
    check_directory_exists "$outdir_dep"
    dep_file="${outdir_dep}/${enlargement}_${cell_line}_AID_Dependent_hotspots_${threshold_num}_${min_count}.bed"

    outdir_indep="${base_dir}/Independent_hotspots"
    check_directory_exists "$outdir_indep"
    indep_file="${outdir_indep}/${enlargement}_${cell_line}_AID_Independent_hotspots_${threshold_num}_${min_count}.bed"
}

# Function to prepare output files by copying header
prepare_outputs() {
    [ -f "$input_file" ] || { echo "Input file not found: $input_file"; exit 1; }

    # Copy the header into the output files
    head -n 1 "$input_file" > "$dep_file"
    head -n 1 "$input_file" > "$indep_file"
}

process_file() {
    echo -e "\nProcessing file: $input_file"
    echo "Dividing by sample: $target_sample"

    # Estrai l'indice della colonna target
    sample_index=$(awk -v sample="$target_sample" '
        BEGIN{FS="\t"}
        NR==1{
            for(i=1;i<=NF;i++){
                if($i==sample){
                    print i;
                    exit;
                }
            }
            exit 1;
        }
    ' "$input_file")

    if [ $? -ne 0 ] || [ -z "$sample_index" ]; then
        echo "Error: Sample $target_sample not found in header!"
        exit 1
    fi

    echo -e "Found sample '$target_sample' at column index: $sample_index\n"

    # Usa awk per separare Dependent/Independent senza loop lento
    awk -v FS="\t" -v OFS="\t" \
        -v idx="$sample_index" \
        -v dep="$dep_file" \
        -v indep="$indep_file" \
        'NR>1 {
            if($idx=="-")
                print >> dep;
            else if($idx=="+")
                print >> indep;
        }' "$input_file"
}


#######################
#        MAIN         #
#######################

start_timer "DepIndep_hotspots_split"
check_args "$@"
build_paths "$@"
prepare_outputs
process_file
check_output_generation "$dep_file"
check_output_generation "$indep_file"
end_timer "DepIndep_hotspots_split"

echo -e "\nFiltering completed!"

