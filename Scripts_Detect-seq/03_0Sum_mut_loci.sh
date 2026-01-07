#!/bin/bash
 
# nohup ./03_0Sum_mut_loci.sh 3 > nohup_03Sum_mut_loc3.out 2>&1 &
# nohup ./03_0Sum_mut_loci.sh 4 > nohup_03Sum_mut_loc4.out 2>&1 &  

###########################################
#              FUNCTIONS                  #
###########################################

source /scratch/scripts/General_Functions.sh

# Funzione che processa ogni file .pmat: conta loci e somma mutazioni
sum_mutations_loci() {
    local file_path="$1"
    local output_locus="$2"
    local output_mut="$3"

    check_file_exists "$file_path"

        local filename=$(basename "$file_path")
        
        # Elaborazione con awk in un solo passaggio
        awk -F '\t' -v fname="$filename" '
        {
            sum += $8
            count += 1
        }
        END {
            print fname "\t" count >> "'"$output_locus"'"
            print fname "\t" sum   >> "'"$output_mut"'"
        }' "$file_path"

        echo -e "\n[INFO] Processed $filename: loci and mutation sum written."
}

###########################################
#                 MAIN                    #
###########################################

# Check argomento soglia iniziale per grafici
if [ $# -ne 1 ]; then
  echo "Usage: $0 <start_threshold>"
  exit 1
fi

start_threshold="$1"

# Directory principale
pmat_dir="/scratch/pmat/Filtered_pmat"

# Record the start time of the entire analysis
start_timer "Sum_mutations_loci_analysis"
echo -e "\n"

# Verifica esistenza
check_directory_exists "$pmat_dir"

# Loop su ciascun sample
for sample_dir in "$pmat_dir"/*; do
    if [ -d "$sample_dir" ]; then
        echo -e "\n-----------------------------------------------------------------\n"
        echo -e "[INFO] Enter the directory: $sample_dir\n"

        filebase=$(basename "$sample_dir")
        echo -e "[INFO] Analysis started on sample: $filebase"

        # Record the start time
        start_timer "$filebase"

        echo -e "\n-------------------------"

        # Funzione interna per filtered_CT e filtered_GA
        for type in "CT" "GA"; do
        echo -e "\n[INFO] Processing type: $type"
            data_type="filtered_${type}"
            data_dir="${sample_dir}/${data_type}"

            if [ -d "$data_dir" ]; then

                # Clean files if they already exist
                for f in \
                    "${data_dir}/Sum_${type}_locnum_${filebase}.txt" \
                    "${data_dir}/Sum_${type}_mutnum_${filebase}.txt"
                do
                    remove_file_if_exists "$f"
                done

                # Create output files
                output_locus="${data_dir}/Sum_${type}_locnum_${filebase}.txt"
                output_mut="${data_dir}/Sum_${type}_mutnum_${filebase}.txt"

                threshold_dirs=($(find "$data_dir" -maxdepth 1 -type d -name "threshold_*" | sort -V))

                for thr_dir in "${threshold_dirs[@]}"; do
                    echo -e "\n[INFO] Processing directory: $thr_dir"
                    for pmat_file in "$thr_dir"/*.pmat; do
                        sum_mutations_loci "$pmat_file" "$output_locus" "$output_mut"
                    done
                done
            else
                echo -e "[ERROR] Directory $data_dir does not exist\n"
            fi
            echo -e "\n"
            check_output_generation "$output_locus"
            check_output_generation "$output_mut"

            # Create directory for graphs if it doesn't exist
            graph_dir="${data_dir}/Graphs"
            mkdir -p "$graph_dir"

            echo -e "\n[INFO] Generating graphs in: $graph_dir"
            # Generate graph for number of loci
            Rscript /scratch/scripts/03_1Sum_mut_loci_graph.R "$output_locus" "$graph_dir" "$start_threshold" "$filebase"

            # Generate graph for number of mutations
            Rscript /scratch/scripts/03_1Sum_mut_loci_graph.R "$output_mut" "$graph_dir" "$start_threshold" "$filebase"

            check_output_generation "$graph_dir/Sum_${type}_locnum_${filebase}_from_${start_threshold}.png"
            check_output_generation "$graph_dir/Sum_${type}_mutnum_${filebase}_from_${start_threshold}.png"

            echo -e "\n-------------------------"

        done
        
    fi
    
    # Record the end time
    echo -e "\n"
    end_timer "$filebase"
    
done

# Record the end time of the entire analysis
end_timer "Sum_mutations_loci_analysis"

echo -e "[INFO] All samples processed successfully.\n"







