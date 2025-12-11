#!/bin/bash

# Docker image name
IMAGE_NAME="detectseqpipe:4"

# Base directory for the project
BASE_DIR="/30tb/home/nuvobea/pmat_and_mpmat/CH12F3"

# Cell line to process
CELL_LINE=("CH12F3")
THRESHOLD=3
MIN_COUNT=5
ENLARGEMENT=0

# Parameters array
PARAMS=("$THRESHOLD" "$MIN_COUNT" "$ENLARGEMENT")

# Samples'order
ORDERED_SAMPLES=("AID_KO_CIT" "UNG_DKO_NS" "UNG_DKO_CIT" "UNG_DKO_CIT_Duv" 
                "UM_TKO_CIT" "UM_TKO_CIT_Duv" "UM_TKO_CIT_Duv_Plain")

# Save original arguments
ORIGINAL_ARGS=("$@")

for cell in "${CELL_LINE[@]}"; do
  echo "Processing cell line: $cell"

  # Restore original arguments at the start of each iteration
  set -- "${ORIGINAL_ARGS[@]}"


    # Host paths
    PMAT_DIR="${BASE_DIR}/pmat_and_mpmat"
    SCRIPTS_DIR="${BASE_DIR}/Scripts_DepIndep"
    GENOME_DIR="${BASE_DIR}/Reference_genome"
    REFERENCE_DIR="${BASE_DIR}/Reference_data"
    DEPINDEP_DIR="${BASE_DIR}/DepIndep_dataset"
    [ -d "${DEPINDEP_DIR}" ] || { echo "$DEPINDEP_DIR doesn't exist, creating it..."; mkdir -p "${DEPINDEP_DIR}"; }


    # Container paths
    CONTAINER_PMAT_DIR="/scratch/pmat"
    CONTAINER_SCRIPTS_DIR="/scratch/scripts"
    CONTAINER_GENOME_DIR="/scratch/genome"
    CONTAINER_REFERENCE_DIR="/scratch/reference_data"
    CONTAINER_DEPINDEP_DIR="/scratch/DepIndep_dataset"


    # Prompt for step if not provided
    if [ -z "$1" ]; then
        echo "Which step do you want to run?"
        echo "1) Shared Hotspots Analysis (01_0SharedHotspots.sh)"
        echo "2) Venn/Euler/Upset Plotting Pipeline (02_SharedHotspots_plots.sh)"
        echo "3) AID Dependent/Independent Hotspots Splitting (03_AID_DepIndephotspots_split.sh)"
        echo "4) Filtering Fragile Sites (04_Remove_fragile_regions.R)"
        echo "5) Ranking AID Dependent/Independent Hotspots (04_0Ranking_AID_DepIndep_hotspots.sh)"
        echo "6) On/Off Dependent/Independent Hotspots Finding (05_OnOff_DepIndep_finding.R)"
        echo "7) Unified Top 200 Hotspots Selection (07_top200hotspots_selection.sh)"
       
        read -p "Enter the number (0/1/2/3/4/5/6/7): " STEP
    else
        STEP="$1"
    fi

    # Shift the step argument to process the remaining as parameters
    shift

    # STEP_ARGS contains the remaining arguments after the step
    STEP_ARGS=("$@")

    # Based on selected step
    case "$STEP" in

        1|"Shared_Hotspots")
            echo "[INFO] Running Shared Hotspots Analysis..."
            # USAGE: $0 1
            # Example:
            # nohup ./00_DepIndep_launch.sh 1 > nohup_01SharedHotspots35_merge0.out 2>&1 &

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/01_0SharedHotspots.sh" "$cell" "${PARAMS[@]}" "${ORDERED_SAMPLES[@]}"
        ;;


        2|"Plots")
            echo "[INFO] Running Venn/Euler/Upset Plotting Pipeline..."
            # USAGE: $0 2 [FLAGS...] (e.g., -e -u)
            # Example:
            # nohup ./00_DepIndep_launch.sh 2 > nohup_02SharedHotspots_plots350.out 2>&1 &


            # Docker execution 
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/02_0SharedHotspots_plots.sh" "$cell" "${PARAMS[@]}" "${STEP_ARGS[@]}"
        ;;


        3|"AID_DepIndep_Split")
            echo "[INFO] Running AID Dependent/Independent Hotspots Splitting..."
            # USAGE: $0 3 <REFERENCE_SAMPLE_NAME>
            # Example:
            # nohup ./00_DepIndep_launch.sh 3 AID_KO_CIT > nohup_03AID_DepIndephotspots_split350.out 2>&1 &
            
            # Requirements: <REFERENCE_SAMPLE_NAME>
            REF_SAMPLE_NAME="${STEP_ARGS[0]}"
            if [ -z "$REF_SAMPLE_NAME" ]; then
                echo "[ERROR] Step 3 requires the reference sample name as the first argument."
                echo "[INFO] Usage: $0 3 <REFERENCE_SAMPLE_NAME>"
                echo "[INFO] Example: $0 3 AID_KO_CIT"
                exit 1
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/03_AID_DepIndephotspots_split.sh" "$cell" "${PARAMS[@]}" "$REF_SAMPLE_NAME"
        ;;


        4|"Ranking_AID_DepIndep")
            echo "[INFO] Running AID Dependent/Independent Hotspots Splitting..."
            # USAGE: $0 5 <ANALYSIS_SAMPLE_NAME> <TYPE> <MODE>
            # <ANALYSIS_SAMPLE_NAME> mandatory
            # <TYPE> optional: --ALL (default), --Dep, --Indep (if nothing entered, both will be done)
            # <MODE> optional: --ALL (default), --mut, --dens (if nothing entered, both will be done)
            # Example:
            # nohup ./00_DepIndep_launch.sh 4 UM_TKO_CIT_Duv > nohup_04Ranking_AID_DepIndep350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 4 UM_TKO_CIT_Duv --Dep > nohup_04Ranking_AID_Dep350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 4 UM_TKO_CIT_Duv --Indep > nohup_04Ranking_AID_Indep350.out 2>&1 &
            
            ANALYSIS_SAMPLE_NAME="${STEP_ARGS[0]}"
            if [ -z "$ANALYSIS_SAMPLE_NAME" ]; then
                echo "[ERROR] Step 4 requires the reference sample name as the first argument."
                echo "[INFO] Usage: $0 4 <ANALYSIS_SAMPLE_NAME>"
                echo "[INFO] Example: $0 4 UM_TKO_CIT_Duv"
                exit 1
            fi

            # Flags are the rest of the arguments (can be empty)
            FLAGS=("${STEP_ARGS[@]:1}")

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/04_0Ranking_AID_DepIndep_hotspots.sh" "$cell" "${PARAMS[@]}" "$ANALYSIS_SAMPLE_NAME" "${FLAGS[@]}"
        ;;

        5|"Filtering_FragileSites")
            echo "[INFO] Running Filtering Fragile Sites..."
            # USAGE: $0 5 <ANALYSIS_SAMPLE_NAME> <TYPE> <MODE>
            # <ANALYSIS_SAMPLE_NAME> mandatory
            # <TYPE> optional: --ALL (default), --Dep, --Indep (if nothing entered, both will be done)
            # <MODE> optional: --ALL (default), --mut, --dens (if nothing entered, both will be done)
            # Example:
            # nohup ./00_DepIndep_launch.sh 5 UM_TKO_CIT_Duv > nohup_05Filtering_FragileSites350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 5 UM_TKO_CIT_Duv --Dep > nohup_05Filtering_FragileSites_Dep350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 5 UM_TKO_CIT_Duv --Indep > nohup_05Filtering_FragileSites_Inde350.out 2>&1 &
            
            ANALYSIS_SAMPLE_NAME="${STEP_ARGS[0]}"
            if [ -z "$ANALYSIS_SAMPLE_NAME" ]; then
                echo "[ERROR] Step 5 requires the reference sample name as the first argument."
                echo "[INFO] Usage: $0 5 <ANALYSIS_SAMPLE_NAME>"
                echo "[INFO] Example: $0 5 UM_TKO_CIT_Duv"
                exit 1
            fi

            # Flags are the rest of the arguments (can be empty)
            FLAGS=("${STEP_ARGS[@]:1}")

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                -v "${REFERENCE_DIR}":"${CONTAINER_REFERENCE_DIR}" \
                "${IMAGE_NAME}" \
                Rscript "${CONTAINER_SCRIPTS_DIR}/05_Remove_fragile_sites.R" "$cell" "${PARAMS[@]}" "$ANALYSIS_SAMPLE_NAME" "${FLAGS[@]}"
        ;;

        6|"OnOff_DepIndep_finding")
            echo "[INFO] Running AID Dependent/Independent Hotspots Splitting..."
            # USAGE: $0 6 <ANALYSIS_SAMPLE_NAME> <TYPE> <MODE>
            # <ANALYSIS_SAMPLE_NAME> mandatory
            # <TYPE> optional: --ALL (default), --Dep, --Indep (if nothing entered, both will be done)
            # <MODE> optional: --ALL (default), --mut, --dens (if nothing entered, both will be done)
            # Example:
            # nohup ./00_DepIndep_launch.sh 6 UM_TKO_CIT_Duv > nohup_06OnOff_Finding_DepIndep350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 6 UM_TKO_CIT_Duv --Dep > nohup_06OnOff_Finding_Dep350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 6 UM_TKO_CIT_Duv --Indep > nohup_06OnOff_Finding_Indep350.out 2>&1 &
            
            ANALYSIS_SAMPLE_NAME="${STEP_ARGS[0]}"
            if [ -z "$ANALYSIS_SAMPLE_NAME" ]; then
                echo "[ERROR] Step 6 requires the reference sample name as the first argument."
                echo "[INFO] Usage: $0 6 <ANALYSIS_SAMPLE_NAME>"
                echo "[INFO] Example: $0 6 UM_TKO_CIT_Duv"
                exit 1
            fi

            # Flags are the rest of the arguments (can be empty)
            FLAGS=("${STEP_ARGS[@]:1}")

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${REFERENCE_DIR}":"${CONTAINER_REFERENCE_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                Rscript "${CONTAINER_SCRIPTS_DIR}/06_OnOff_DepIndep_finding.R" "$cell" "${PARAMS[@]}"  "$ANALYSIS_SAMPLE_NAME" "${FLAGS[@]}"
        ;;


        7|"Top200_Selection")
            echo "[INFO] Running Unified Top 200 Hotspots Selection..."
            # USAGE: $0 7 <ANALYSIS_SAMPLE_NAME> <TYPE> <MODE>
            # <ANALYSIS_SAMPLE_NAME> mandatory
            # <TYPE> optional: --ALL (default), --Dep, --Indep
            # <MODE> optional: --ALL (default), --mut, --dens
            # Example:
            # nohup ./00_DepIndep_launch.sh 7 UM_TKO_CIT_Duv > nohup_07top200_selection350.out 2>&1 &
 
            ANALYSIS_SAMPLE_NAME="${STEP_ARGS[0]}"
            if [ -z "$ANALYSIS_SAMPLE_NAME" ]; then
             echo "[ERROR] Step 7 requires the analysis sample name as the first argument."
             echo "[INFO] Usage: $0 7 <ANALYSIS_SAMPLE_NAME>"
             echo "[INFO] Example: $0 7 UM_TKO_CIT_Duv_CLEAN"
             exit 1
            fi

            # Flags are the rest of the arguments (can be empty)
            FLAGS=("${STEP_ARGS[@]:1}")

            # Docker execution
            docker run --rm -i \
             --user $(id -u):$(id -g) \
             -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
             -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
             -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
            "${IMAGE_NAME}" \
            bash "${CONTAINER_SCRIPTS_DIR}/07_top200hotspots_selection.sh" "$cell" "${PARAMS[@]}" "$ANALYSIS_SAMPLE_NAME" "${FLAGS[@]}"
        ;;

        *)
            echo "Invalid option: $STEP"
            echo "Usage: $0 [0|1|2|3|4|5|6|7] [args...]"
        ;;
    esac
done