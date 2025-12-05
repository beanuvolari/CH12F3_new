#!/bin/bash

# nohup ./33_4Ranking_Dep_Indep_hotspots200.sh 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN CH12F3 > nohup_ranking200_Dep_UNG_KO_blackfilt.out 2>&1 &
# nohup ./33_4Ranking_Dep_Indep_hotspots200.sh 3 5 0 Indep UM_TKO_CIT_Duv_CLEAN CH12F3 > nohup_ranking200_Indep_UNG_KO.out 2>&1 &

############################
#        FUNCTIONS         #
############################

source ./General_Functions.sh

# Check input parameters
function check_params() {
    if [ "$#" -lt 5 ]; then
        echo "[ERROR] You must provide at least 5 parameters:"
        echo "[INFO] Usage: $0 threshold_num min_count enlargement type sample"
        echo "[INFO] type should be 'Dep' or 'Indep'"
        exit 1
    fi
}

# Extract top 200 hotspots from the given input file and write to output
function process_hotspots_file() {
    local input_file=$1
    local output_file=$2

    if [ ! -f "$input_file" ]; then
        echo "[ERROR] Input file not found: $input_file"
        exit 1
    fi

    {
        head -n 1 "$input_file"      # Header line
        head -n 201 "$input_file" | tail -n +2  # First 200 data rows (excluding header)
    } > "$output_file"

    echo "[INFO] Output file saved: $output_file"
}


####################
#       MAIN       #
####################


# Main workflow function
function main() {
    check_params "$@"

    local threshold_num="$1"
    local min_count="$2"
    local enlargement="$3"
    local type="$4"
    local sample="$5"
    local cell_line="$6"

    local dens_muts=("density" "nmutations")
    start_timer "Ranking_Dep_Indep_hotspots_script"

    echo "[INFO] Processing started at: $start_time_human"
    echo

    for dens_mut in "${dens_muts[@]}"; do
        local base_dir="./Dep_Indep_dataset/Threshold_${threshold_num}/Enlargement_${enlargement}/${sample}/Ranking_${type}_hotspots/${dens_mut}"

        if [ ! -d "$base_dir" ]; then
            echo "[ERROR] Directory not found: $base_dir"
            exit 1
        fi

        local input_file="${base_dir}/${enlargement}_${sample}_${type}_${cell_line}_Ranked_hotspots_${dens_mut}_${threshold_num}${min_count}_blackfilt.tsv"
        local output_file="${base_dir}/${enlargement}_${sample}_${type}_${cell_line}_Ranked200hotspots_${dens_mut}_${threshold_num}${min_count}_blackfilt.tsv"

        echo "[INFO] Processing file: $input_file"
        process_hotspots_file "$input_file" "$output_file"
        echo
    done

    end_timer "Ranking_Dep_Indep_hotspots_script"
    echo "[INFO] Process completed"
}

# Start script execution
main "$@"
