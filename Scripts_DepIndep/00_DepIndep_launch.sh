#!/bin/bash

# Docker image name
IMAGE_NAME="detectseqpipe:4"

# Base directory for the project
BASE_DIR="/30tb/home/nuvobea/pmat_and_mpmat/CH12F3"

# Cell line to process
CELL_LINE=("CH12F3")

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
        echo "4) Ranking AID Dependent/Independent Hotspots (04_0Ranking_AID_DepIndep_hotspots.sh)"
       
        read -p "Enter the number (0/1/2/3/4): " STEP
    else
        STEP="$1"
    fi

    # Shift the step argument to process the remaining as parameters
    shift

    # Based on selected step
    case "$STEP" in

        1|"Shared_Hotspots")
            echo "[INFO] Running Shared Hotspots Analysis..."

            # nohup ./00_DepIndep_launch.sh 1 3 5 0 > nohup_01SharedHotspots35_merge0.out 2>&1 &

            # Parameters:
            if [ $# -lt 3 ]; then
                echo "[INFO] You must specify threshold, min_count and enlargement."
                read -p "Enter thresholds (e.g. 3): " THRESHOLDS
                read -p "Enter min_count (e.g. 5): " MIN_COUNT
                read -p "Enter enlargement (e.g.0, 1, 2 or 3): " ENLARGEMENT
                ARGS=("$THRESHOLDS" "$MIN_COUNT" "$ENLARGEMENT")
            else
                ARGS=("$@")
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/01_0SharedHotspots.sh" "$cell" "${ARGS[@]}" "${ORDERED_SAMPLES[@]}"
        ;;


        2|"Plots")
            echo "[INFO] Running Venn/Euler/Upset Plotting Pipeline..."

            # nohup ./00_DepIndep_launch.sh 2 3 5 0 > nohup_02SharedHotspots_plots350.out 2>&1 &
            
            # Logic for arguments specific to the plotting wrapper
            # It expects: <threshold> [OPTIONS like -e, -v, ALL]
            if [ $# -lt 3 ]; then
                echo "[INFO] You must specify threshold, min_count and enlargement."
                read -p "Enter thresholds (e.g. 3): " THRESHOLDS
                read -p "Enter min_count (e.g. 5): " MIN_COUNT
                read -p "Enter enlargement (e.g.0, 1, 2 or 3): " ENLARGEMENT
                read -p "Enter flags (optional, e.g. 'ALL' or '-e -u'. Press Enter for default ALL): " FLAGS
                # Combine them into an array
                ARGS=("$THRESHOLDS" "$MIN_COUNT" "$ENLARGEMENT" "$FLAGS")
            else
                # If passed via command line (e.g. ./script.sh 2 10 -e)
                ARGS=("$@")
            fi

            # Docker execution 
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/02_0SharedHotspots_plots.sh" "$cell" "${ARGS[@]}"
        ;;


        3|"AID_DepIndep_Split")
            echo "[INFO] Running AID Dependent/Independent Hotspots Splitting..."

            # nohup ./00_DepIndep_launch.sh 3 3 5 0 AID_KO_CIT > nohup_03AID_DepIndephotspots_split350.out 2>&1 &
            
            # Parameters:
            if [ $# -lt 3 ]; then
                echo "[INFO] You must specify threshold, min_count, enlargement and the sample needed for the analysis."
                read -p "Enter thresholds (e.g. 3): " THRESHOLDS
                read -p "Enter min_count (e.g. 5): " MIN_COUNT
                read -p "Enter enlargement (e.g.0, 1, 2 or 3): " ENLARGEMENT
                read -p "Enter the target sample name (e.g. AID_KO_CIT): " SAMPLE_NAME
                # Combine them into an array
                ARGS=("$THRESHOLDS" "$MIN_COUNT" "$ENLARGEMENT" "$SAMPLE_NAME")
            else
                ARGS=("$@")
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/03_AID_DepIndephotspots_split.sh" "$cell" "${ARGS[@]}"
        ;;

        4|"Ranking_AID_DepIndep")
            echo "[INFO] Running AID Dependent/Independent Hotspots Splitting..."

            # nohup ./00_DepIndep_launch.sh 4 3 5 0 UM_TKO_CIT_Duv --Dep > nohup_04Ranking_AID_Dep350.out 2>&1 &
            # nohup ./00_DepIndep_launch.sh 4 3 5 0 UM_TKO_CIT_Duv --Indep > nohup_04Ranking_AID_Indep350.out 2>&1 &
            
            # Parameters:
            if [ $# -lt 3 ]; then
                echo "[INFO] You must specify threshold, min_count, enlargement, the sample, the type and the mode (not mandatory)."
                read -p "Enter thresholds (e.g. 3): " THRESHOLDS
                read -p "Enter min_count (e.g. 5): " MIN_COUNT
                read -p "Enter enlargement (e.g.0, 1, 2 or 3): " ENLARGEMENT
                read -p "Enter the target sample name (e.g. AID_KO_CIT): " SAMPLE_NAME
                read -p "Enter the type (e.g. --Dep or --Indep): " TYPE
                read -p "Enter the mode [--mut|--dens] (if you enter nothing both will be done): " MODE
                # Combine them into an array
                ARGS=("$THRESHOLDS" "$MIN_COUNT" "$ENLARGEMENT" "$SAMPLE_NAME" "$TYPE" "$MODE")
            else
                ARGS=("$@")
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PMAT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${DEPINDEP_DIR}":"${CONTAINER_DEPINDEP_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/04_0Ranking_AID_DepIndep_hotspots.sh" "$cell" "${ARGS[@]}"
        ;;

        *)
            echo "Invalid option: $STEP"
            echo "Usage: $0 [0|1|2|3|4] [args...]"
        ;;
    esac
done