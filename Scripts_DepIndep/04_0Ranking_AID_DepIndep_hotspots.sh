#!/bin/bash
# Example usage:
# nohup ./04_0Ranking_AID_DepIndep_hotspots.sh CH12F3 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN Dep > nohup_ranking_Dep_UM_TKO_CIT_Duv.out 2>&1 &
# nohup ./04_0Ranking_AID_DepIndep_hotspots.sh CH12F3 3 5 0 Indep UM_TKO_CIT_Duv_CLEAN Indep > nohup_ranking_Indep_UM_TKO_CIT_Duv.out 2>&1 &

############################
#        FUNCTIONS         #
############################

source /scratch/scripts/General_Functions.sh

check_args() {
    if [ "$#" -lt 6 ]; then
        echo "Error: Missing arguments."
        echo "Usage: $0 <cell_line> <threshold_num> <min_count> <enlargement> <sample>  [ALL|--Dep|--Indep] [ALL|--mut|--dens]"
        exit 1
    fi
}



build_paths() {
    cell_line="$1"
    threshold_num="$2"
    min_count="$3"
    enlargement="$4"
    sample="$5"            # es. UM_TKO_CIT_Duv_CLEAN
    type="$6"


    base_dir="/scratch/DepIndep_dataset/Threshold_${threshold_num}/Enlargement_${enlargement}"

    # Hotspots dataset (da step precedente, ora con cell_line variabile)
    hotspots_file="${base_dir}/${enlargement}_${cell_line}_AID_${type}endent_hotspots_${threshold_num}${min_count}.bed"
    [ -f "$hotspots_file" ] || { echo "Hotspots file not found: $hotspots_file"; exit 1; }
    echo "Hotspots file: $hotspots_file"
    
    # Mutations dataset
    mutations_dir="/scratch/pmat/Filtered_pmat"
    mutations_file="${mutations_dir}/${sample}/merged_CTGA/threshold_${threshold_num}/${sample}_merged_CTGA_${threshold_num}_sorted.bed"
    [ -f "$mutations_file" ] || { echo "Mutations file not found: $mutations_file"; exit 1; }
    echo "Mutations file: $mutations_file"

    # Output path
    output_path="${base_dir}/${sample}/Ranking_${type}_hotspots"
    check_directory_exists "$output_path"

    # Output files (created only if enabled)
    if [[ "$RUN_DENS" == true ]]; then
        check_directory_exists "$output_path/density"
        output_density="${output_path}/density/${enlargement}_${sample}_${type}_${cell_line}_Ranked_hotspots_density_${threshold_num}${min_count}.tsv"
    else
        output_density=""
    fi

    if [[ "$RUN_MUT" == true ]]; then
        check_directory_exists "$output_path/nmutations"
        output_mutation="${output_path}/nmutations/${enlargement}_${sample}_${type}_${cell_line}_Ranked_hotspots_nmutations_${threshold_num}${min_count}.tsv"
    else
        output_mutation=""
    fi

}

run_ranking() {
    echo -e "\nRunning Ranking:"
    echo "Hotspots:   $hotspots_file"
    echo "Mutations:  $mutations_file"
    echo "Sample:     $sample"
    echo "Cell line:  $cell_line"
    echo "Type:       $type"
    echo "Params:     enlargement=$enlargement, threshold=$threshold_num, min_count=$min_count"

    Rscript ./04_1Ranking_AID_DepIndep_hotspots.R \
        "$sample" "$hotspots_file" "$mutations_file" \
        "$output_density" "$output_mutation"

    echo -e "\nGenerated outputs:"
    [[ -n "$output_density" ]] && echo " - $output_density"
    [[ -n "$output_mutation" ]] && echo " - $output_mutation"
}


############################
#   DEFAULT FLAG VALUES    #
############################

# Default: execute everything unless the user restricts it
RUN_DEP=true
RUN_INDEP=true
RUN_DENS=true
RUN_MUT=true

############################
#    PARSE OPTIONAL FLAGS  #
############################
parse_flags() {

    for arg in "$@"; do
        case "$arg" in
            --Dep)
                RUN_DEP=true
                RUN_INDEP=false
                ;;
            --Indep)
                RUN_DEP=false
                RUN_INDEP=true
                ;;
            --mut)
                RUN_MUT=true
                RUN_DENS=false
                ;;
            --dens)
                RUN_MUT=false
                RUN_DENS=true
                ;;
        esac
    done

    echo "[INFO] Flags:"
    echo "       RUN_DEP = $RUN_DEP"
    echo "       RUN_INDEP = $RUN_INDEP"
    echo "       RUN_DENS = $RUN_DENS"
    echo "       RUN_MUT = $RUN_MUT"
}



####################
#       MAIN       #
####################

check_args "$@"

# Extract required positional params
cell_line="$1"
threshold_num="$2"
min_count="$3"
enlargement="$4"
sample="$5"


# Parse optional flags
parse_flags "$@"


start_timer "Ranking_Dep/Indep_hotspots"

# Loop through Dep and Indep depending on flags
for type in Dep Indep; do

    [[ "$type" == "Dep" && "$RUN_DEP" == false ]] && continue
    [[ "$type" == "Indep" && "$RUN_INDEP" == false ]] && continue

    build_paths "$cell_line" "$threshold_num" "$min_count" "$enlargement" "$sample" "$type"
    run_ranking
done

end_timer "Ranking_Dep/Indep_hotspots"

echo -e "\nProcess completed successfully.\n"
