#!/usr/bin/env bash

# nohup ./33_5Takingfirst200hotspots_filtered.sh 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN density CH12F3 > nohup_33_5Takingfirst200hotspotsfilt_Depdens.out 2>&1 &
# nohup ./33_5Takingfirst200hotspots_filtered.sh 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_33_5Takingfirst200hotspotsfilt_Depnmut.out 2>&1 &

##########################
#       FUNCTIONS        #
##########################


filter_rank() {
    local input="$1"
    local output="$2"
    local threshold="${3:-200}"   # default 200 se non specificato

    # Trova la colonna "Rank" nell'intestazione
    local col
    col=$(head -1 "$input" | tr '\t' '\n' | grep -n '^Rank$' | cut -d: -f1)

    if [ -z "$col" ]; then
        echo "[ERROR]: 'Rank' column not found in $input"
        return 2
    fi

    # Filtra le righe (mantiene intestazione + Rank <= soglia)
    awk -v c="$col" -v thr="$threshold" -F'\t' 'NR==1 || $c <= thr' "$input" > "$output"

    echo "[INFO] Filtered file saved in $output (soglia = $threshold)"
}


##########################
#          MAIN          #
##########################

threshold_num="$1"
min_count="$2"
enlargement="$3"
type="$4"              # Dep / Indep
sample="$5"            # es. Primary_B_cells_UNG_KO_CLEAN
dens_mut="$6"
cell_line="$7"         # es. PrimaryB

base_dir="./Dep_Indep_dataset/Threshold_${threshold_num}/Enlargement_${enlargement}/${sample}/Ranking_${type}_hotspots/${dens_mut}"



input_filt="${base_dir}/${enlargement}_${sample}_${type}_${cell_line}_Ranked_hotspots_${dens_mut}_${threshold_num}${min_count}_blackfilt.tsv"
input_black="${base_dir}/${enlargement}_${sample}_${type}_${cell_line}_Deleted_hotspots_${dens_mut}_${threshold_num}${min_count}_blacklist.tsv"

output_dir="${base_dir}/Filtered_and_blacklist_top200"
mkdir -p "$output_dir"
output_filt="${output_dir}/${enlargement}_${sample}_${type}_${cell_line}_Ranked_hotspots_${dens_mut}_${threshold_num}${min_count}_200blackfilt.tsv"
output_black="${output_dir}/${enlargement}_${sample}_${type}_${cell_line}_Deleted_hotspots_${dens_mut}_${threshold_num}${min_count}_200blacklist.tsv"

filter_rank "$input_filt" "$output_filt" 200
filter_rank "$input_black" "$output_black" 200

echo "[INFO] All filtered files have been saved in $output_dir"
