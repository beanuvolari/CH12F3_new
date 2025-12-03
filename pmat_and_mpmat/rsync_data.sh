#!/bin/bash
# nohup ./rsync_data.sh > nohup_rsync.out 2>&1 &
# This script is used to synchronize data files from a source directory to a destination directory.

array1=("UM_TKO_CIT_CLEAN.pmat" "UM_TKO_CIT_Duv_CLEAN.pmat" "UM_TKO_CIT_Duv_Plain_CLEAN.pmat")
array2=("UNG_DKO_NS_CLEAN.pmat" "UNG_DKO_CIT_CLEAN.pmat" "UNG_DKO_CIT_Duv_CLEAN.pmat")

log_time() {
    local start=$1
    local end=$2
    local file=$3
    local duration=$((end - start))
    echo "Finished copying $file in $duration seconds."
}

# Copy files TKO
echo -e "Starting rsync for TKO files..."
for file1 in "${array1[@]}"; do
    echo "Processing file: $file1"
    start_time=$(date +%s)
    rsync -av /30tb/alessandri/detectSeq2_V1/Exp1/pmat_and_mpmat/"$file1" \
        /30tb/home/nuvobea/pmat_and_mpmat/CH12F3/Original_pmat/ \
   
    end_time=$(date +%s)
    log_time "$start_time" "$end_time" "$file1"
done

# Copy files DKO
echo -e "Starting rsync for DKO files..."
for file2 in "${array2[@]}"; do
    echo "Processing file: $file2"
    start_time=$(date +%s)
    rsync -av /30tb/alessandri/detectSeq2_V1/Exp1/pmat_and_mpmat/"$file2" \
        /30tb/home/nuvobea/pmat_and_mpmat/CH12F3/Original_pmat/ \
   
    end_time=$(date +%s)
    log_time "$start_time" "$end_time" "$file2"
done

# Ultimo file singolo
file3="AID_KO_CIT_merged_sorted_CLEAN.pmat"
echo "Processing file: $file3"
start_time=$(date +%s)
rsync -av /30tb/alessandri/detectSeq2/Exp1/pmat_and_mpmat/"$file3" \
    /30tb/home/nuvobea/pmat_and_mpmat/CH12F3/Original_pmat/ \
  
end_time=$(date +%s)
log_time "$start_time" "$end_time" "$file3"
echo "All files have been processed successfully."
echo "Script execution completed."


# nohup rsync -av --progress /30tb/alessandri/detectSeq2_V1/Exp1/pmat_and_mpmat/UM_TKO_CIT_CLEAN.pmat /30tb/home/nuvobea/pmat_and_mpmat/CH12F3/Original_pmat/ > rsync.log 2>&1 &