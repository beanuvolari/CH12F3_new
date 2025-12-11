#!/usr/bin/env bash

# This script performs two main steps for all specified Type/Mode combinations:
# 1. Positional extraction (head/tail) of the Top 200 hotspots from the FULL and FILTERED files (Script 07 logic).
# 2. Content-based extraction (awk) of hotspots with Rank <= 200 from the FILTERED and BLACKLIST files (Script 08 logic).

# Example Usages:
# nohup ./07_top200hotspots_selection.sh CH12F3 3 5 0 UM_TKO_CIT_Duv > nohup_07top200hotspots_selection.out 2>&1 &

############################
#       CONFIGURATION      #
############################

source /scratch/scripts/General_Functions.sh

# --- DEFAULT FLAG VALUES ---
RUN_DEP=true
RUN_INDEP=true
RUN_DENS=true
RUN_MUT=true

############################
#         FUNCTIONS        #
############################

# Check required positional parameters
check_params() {
    if [ "$#" -lt 5 ]; then
        echo "[ERROR] You must provide at least 5 positional parameters."
        echo "[INFO] Usage: $0 cell_line threshold_num min_count enlargement sample [optional flags...]"
        exit 1
    fi
}

# --- PARSE OPTIONAL FLAGS ---
parse_flags() {
    # Initialization and logic remain the same as defined previously
    RUN_DEP=false; RUN_INDEP=false; RUN_MUT=false; RUN_DENS=false
    TYPE_SET=false; MODE_SET=false; ALL_SEEN=false
    
    for arg in "$@"; do
        case "$arg" in
            --ALL) ALL_SEEN=true ;;
            --Dep) RUN_DEP=true; RUN_INDEP=false; TYPE_SET=true ;;
            --Indep) RUN_DEP=false; RUN_INDEP=true; TYPE_SET=true ;;
            --mut) RUN_MUT=true; RUN_DENS=false; MODE_SET=true ;;
            --dens) RUN_MUT=false; RUN_DENS=true; MODE_SET=true ;;
        esac
    done

    if [[ "$TYPE_SET" == false ]]; then RUN_DEP=true; RUN_INDEP=true; fi
    if [[ "$MODE_SET" == false ]]; then RUN_MUT=true; RUN_DENS=true; fi

    echo "[INFO] Flags set for execution:"
    echo "[INFO] RUN_DEP=$RUN_DEP  RUN_INDEP=$RUN_INDEP  RUN_MUT=$RUN_MUT  RUN_DENS=$RUN_DENS"
}

# --- SCRIPT 07 LOGIC: Positional Extraction ---
# Extracts top 200 hotspots (positional selection)
process_top_200_hotspots() {
    local input_file="$1"
    local output_file="$2"
    local run_type="$3"
    local run_mode="$4"
    
    if [ ! -f "$input_file" ]; then
        echo "[WARNING] Input file not found for ${run_type}/${run_mode} (Positional): $input_file"
        return 1
    fi

    echo "[INFO] Starting positional extraction (rows 1-201) for ${run_type}/${run_mode}"
    {
        head -n 1 "$input_file"     # Header line
        head -n 201 "$input_file" | tail -n +2  # First 200 data rows (excluding header)
    } > "$output_file"

    echo "[INFO] Positional Output saved: $output_file"
    return 0
}

# --- SCRIPT 08 LOGIC: Content-Based Filtering ---
# Filters rows where the value in the 'Rank' column is <= threshold
filter_rank() {
    local input="$1"
    local output="$2"
    local run_type="$3"
    local run_mode="$4"
    local threshold="${5:-200}"    # Default 200 if not specified

    if [ ! -f "$input" ]; then
        echo "[WARNING] Input file not found for ${run_type}/${run_mode} (Rank Filter): $input"
        return 1
    fi

    # Find the column number of "Rank" in the header
    local col
    col=$(head -1 "$input" | tr '\t' '\n' | grep -n '^Rank$' | cut -d: -f1)

    if [ -z "$col" ]; then
        echo "[ERROR]: 'Rank' column not found in $input. Skipping Rank filtering."
        return 2
    fi

    # Filter rows (keep header + Rank <= threshold)
    awk -v c="$col" -v thr="$threshold" -F'\t' 'NR==1 || $c <= thr' "$input" > "$output"

    echo "[INFO] Rank Filtered file saved in $output (Threshold = $threshold)"
    return 0
}


####################
#        MAIN      #
####################

check_params "$@"

cell_line="$1"; threshold_num="$2"; min_count="$3"; enlargement="$4"; sample="$5"

parse_flags "$@"

start_timer "Complete_Unified_Hotspot_Processing_script"
echo "[INFO] Processing started at: $start_time_human"
echo

types=()
if [ "$RUN_DEP" == true ]; then types+=("Dep"); fi
if [ "$RUN_INDEP" == true ]; then types+=("Indep"); fi

modes=()
if [ "$RUN_MUT" == true ]; then modes+=("nmutations"); fi
if [ "$RUN_DENS" == true ]; then modes+=("density"); fi


# Double loop over selected Types and Modes
for current_type in "${types[@]}"; do
    for current_mode in "${modes[@]}"; do
        
        echo "------------------------------------------------------------------"
        echo "[INFO] Starting processing for TYPE: ${current_type} and MODE: ${current_mode}"

        # --- Directory Setup ---
        base_dir="/scratch/DepIndep_dataset/Threshold_${threshold_num}/Enlargement_${enlargement}/${current_type}endent_hotspots/${sample}/Ranking_${current_type}_hotspots/${current_mode}"
        filter_dir="${base_dir}/Filtered_hotspots"
        output_awk_dir="${filter_dir}/Filtered_and_blacklist_top200"

        if [ ! -d "$base_dir" ]; then
            echo "[ERROR] Base directory not found: $base_dir. Skipping."
            continue
        fi
        
        # Ensure the filter directory exists for writing files (Script 07 outputs)
        mkdir -p "$filter_dir"
        # Ensure the final output directory exists (Script 08 outputs)
        mkdir -p "$output_awk_dir"


        # --- FILES: 1. FULL Set (Script 07 Positional Logic) ---
        input_full="${base_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Ranked_hotspots_${current_mode}_${threshold_num}${min_count}.tsv"
        output_full_pos="${base_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Ranked200hotspots_${current_mode}_${threshold_num}${min_count}.tsv"
        
        process_top_200_hotspots "$input_full" "$output_full_pos" "$current_type" "$current_mode"


        # --- FILES: 2. FILTERED Set (Script 07 Positional Logic & Script 08 Rank Logic) ---
        input_filt="${filter_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Ranked_hotspots_${current_mode}_${threshold_num}${min_count}_blackfilt.tsv"
        output_filt_pos="${filter_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Ranked200hotspots_${current_mode}_${threshold_num}${min_count}_blackfilt.tsv"
        output_filt_awk="${output_awk_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Ranked_hotspots_${current_mode}_${threshold_num}${min_count}_200blackfilt.tsv"

        # Script 07 Positional Logic on Filtered file
        process_top_200_hotspots "$input_filt" "$output_filt_pos" "$current_type" "$current_mode"
        
        # Script 08 Rank Logic on Filtered file
        filter_rank "$input_filt" "$output_filt_awk" "$current_type" "$current_mode" 200


        # --- FILES: 3. BLACKLIST Set (Script 08 Rank Logic) ---
        input_black="${filter_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Deleted_hotspots_${current_mode}_${threshold_num}${min_count}_blacklist.tsv"
        output_black_awk="${output_awk_dir}/${enlargement}_${sample}_${current_type}_${cell_line}_Deleted_hotspots_${current_mode}_${threshold_num}${min_count}_200blacklist.tsv"

        # Script 08 Rank Logic on Blacklist file
        filter_rank "$input_black" "$output_black_awk" "$current_type" "$current_mode" 200


        echo "[INFO] Combination ${current_type}/${current_mode} finished. All 4 Top 200 files created."

    done
done

end_timer "Complete_Unified_Hotspot_Processing_script"
echo "[INFO] Process completed successfully."