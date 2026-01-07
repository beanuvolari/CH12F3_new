#!/bin/bash

##################
# SET PARAMETERS #
##################

# Docker image name
IMAGE_NAME="detectseqpipe:4"
IMAGE_NAME2="repbioinfo/detectseq"

# Base directory for the project
BASE_DIR="/30tb/home/nuvobea/pmat_and_mpmat/CH12F3"

# Cell line to process
CELL_LINE=("CH12F3")

# Adapters
ADAPT1="AGATCGGAAGAGCACACGT"
ADAPT2="AGATCGGAAGAGCGTCGTG"



##################
# START SCRIPTS  #
##################

# --- HELP ---
show_help() {
    echo "=========================================================================="
    echo " DETECT-SEQ PIPELINE LAUNCHER"
    echo "=========================================================================="
    echo "Usage: ./00_Pipeline_launch.sh [STEP] [ARGUMENTS]"
    echo ""
    echo "STEPS:"
    echo "  1) fastq-to-pmat       Generate PMAT files from FASTQ."
    echo "                         Args: --all/nothing OR --single [sample_R1.fastq.gz ...]"
    echo ""
    echo "  2) filtering           Run PMAT filtering."
    echo "                         Args: --all/nothing OR --part1, --part2 [thresholds]"
    echo ""
    echo "  3) sum-mut-loci        Sum mutation loci and generate graphs."
    echo "                         Args: [start_threshold] (e.g., 3 or 6)"
    echo ""
    echo "  4) merging-ctga        Merge CTGA data."
    echo "                         Args: [thresholds (optional)] (e.g., 3-5)"
    echo ""
    echo "  5) rider-processing    Process Rider data."
    echo "                         Args: [mode: --single|--merge] [min_count] [thresholds (optional)]"
    echo ""
    echo "EXAMPLES:"
    echo "  nohup ./00_Pipeline_launch.sh 1 --all > log.out 2>&1 &"
    echo "  nohup ./00_Pipeline_launch.sh 2 --part2 3-5 > log.out 2>&1 &"
    echo "  nohup ./00_Pipeline_launch.sh 5 --merge 5 3-5 > log.out 2>&1 &"
    echo "=========================================================================="
}

# Check if no arguments are provided or help is requested
if [ -z "$1" ] || [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    show_help
    exit 0
fi
# --- END HELP ---

# --- MAIN PIPELINE ---
# Save original arguments
ORIGINAL_ARGS=("$@")

for cell in "${CELL_LINE[@]}"; do
  echo "Processing cell line: $cell"

  # Restore original arguments at the start of each iteration
  set -- "${ORIGINAL_ARGS[@]}"

  #exp_dir="${BASE_DIR}/Results/${cell}"
  #mkdir -p "${exp_dir}"

    # Host paths
   # FASTQ_DIR="${BASE_DIR}/Original_fastq/${cell}/merged_fastq"
    PROJECT_DIR="${BASE_DIR}/pmat_and_mpmat"
    SCRIPTS_DIR="${BASE_DIR}/Scripts_Detect-seq"
    GENOME_DIR="${BASE_DIR}/Reference_genome"

    # Directories needed for step 4 and 5
    RIDER_DIR="${BASE_DIR}/Rider"
    REFERENCE_DIR="${BASE_DIR}/Reference_data"

    # Container paths
    #CONTAINER_FASTQ_DIR="/scratch/raw_fastq:ro"
    CONTAINER_PMAT_DIR="/scratch/pmat"
    CONTAINER_SCRIPTS_DIR="/scratch/scripts"
    CONTAINER_GENOME_DIR="/scratch/genome"
    CONTAINER_RIDER_DIR="/scratch/rider"
    CONTAINER_REFERENCE_DIR="/scratch/reference_data"


    # Check if no step is provided (Redundant but kept for structure)
    if [ -z "$1" ]; then
        echo "[ERROR] No step provided. Use --help to see usage."
        echo "[INFO] Usage: ./00_Pipeline_launch.sh --help or -h"
        exit 1
    fi

    STEP="$1"

    # Shift the step argument to process the remaining as parameters
    shift

    # Based on selected step
    case "$STEP" in

        1|"fastq-to-pmat")
            echo "[INFO] Launching 01_pmat_generation.sh..."
            # USAGE: $0 1 [FLAGS: --all/nothing | --single sample_R1.fastq.gz ...]
            # Examples:
            # nohup ./00_Pipeline_launch.sh 1 --all > nohup_01_all_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 1 > nohup_01_all_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 1 --single sample_R1.fastq.gz sample2_R1.fastq.gz > nohup_01_single_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 1 --single AID_KO_CIT_merged_R1.fastq.gz > nohup_01_single_AID_KO_CIT_pmat.out 2>&1 &

            # Default behavior: if not specified elaborate all samples in raw_fastq
            MODE="--all"
            SAMPLE_ARGS=()

           # Parse mode
           MODE="${1:---all}" # Default a --all se vuoto
            if [ "$MODE" == "--single" ]; then
                shift
                if [ $# -eq 0 ]; then
                    echo "[ERROR] Step $STEP in --single mode requires at least one sample name."
                    echo "[INFO] Usage: $0 $STEP --single sample_R1.fastq.gz"
                    echo -e "\n[INFO] Showing help:\n"
                    show_help
                    exit 1
                fi
                SAMPLE_ARGS=("$@")
            else
                SAMPLE_ARGS=()
                [ "$MODE" == "--all" ] && shift
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${GENOME_DIR}":"${CONTAINER_GENOME_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                -v "${FASTQ_DIR}":"${CONTAINER_FASTQ_DIR}" \
                -v "${exp_dir}":/"${CONTAINER_EXP_DIR}" \
                "$IMAGE_NAME2" \
                bash "$CONTAINER_SCRIPTS_DIR/01_pmat_generation.sh" "$ADAPT1" "$ADAPT2" "$MODE" "${SAMPLE_ARGS[@]}"
            ;;

        2|"filtering")
            echo "[INFO] Launching 02_Filtering_pmat.sh..."
            # USAGE: $0 2 [FLAGS: --all/nothing | --part1 | --part2] <THRESHOLDS>
            # Examples:

            # nohup ./00_Pipeline_launch.sh 2 3-5 > nohup_02_3-5_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 6-10 > nohup_02_6-10_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 --part2 6 7 8 > nohup_02_678_part2.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 --part2 6-10 > nohup_02_6-10_part2.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 --part1 > nohup_02_part1.out 2>&1 &

           if [ -z "$1" ]; then
                echo "[ERROR] Step $STEP requires thresholds to filter."
                echo "[INFO] Usage: $0 $STEP [FLAGS: (--all, --part1, --part2)] <THRESHOLDS>"
                echo "[INFO] Example: $0 $STEP 3-5  or  $0 $STEP --part2 6 7 8"
                echo -e "\n[INFO] Showing help:\n"
                show_help
                exit 1
            fi
        
            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PROJECT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/02_Filtering_pmat.sh" "$@"
            ;;

        3|"sum"|"sum-mut-loci")
            echo "[INFO] Launching 03_0Sum_mut_loci.sh..."
            # USAGE: $0 3 <THRESHOLD>
            # Examples:
            # nohup ./00_Pipeline_launch.sh 3 6 > nohup_03_Sum_mut_loci_6.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 3 3 > nohup_03_Sum_mut_loci_3.out 2>&1 &
            
            START_THRESHOLD="$1"
            if [ -z "$START_THRESHOLD" ]; then
                echo "[ERROR] Step $STEP requires a start_threshold argument."
                echo "[INFO] Usage: $0 $STEP <THRESHOLD>"
                echo "[INFO] Example: $0 $STEP 3"
                echo -e "\n[INFO] Showing help:\n"
                show_help
                exit 1
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "$PROJECT_DIR":"$CONTAINER_PMAT_DIR" \
                -v "$SCRIPTS_DIR":"$CONTAINER_SCRIPTS_DIR" \
                "$IMAGE_NAME" \
                bash "$CONTAINER_SCRIPTS_DIR/03_0Sum_mut_loci.sh" "$START_THRESHOLD"
            ;;

        4|"merge"|"merging-ctga")
            echo "Launching 04_Merging_CTGA.sh script..."
            # USAGE: $0 4 [THRESHOLDS/nothing]
            # Examples:
            # nohup ./00_Pipeline_launch.sh 4 > nohup_04_Merging_CTGA.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 4 3-5 > nohup_04_Merging_CTGA_3-5.out 2>&1 &


            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "$PROJECT_DIR":"$CONTAINER_PMAT_DIR" \
                -v "$SCRIPTS_DIR":"$CONTAINER_SCRIPTS_DIR" \
                -v "$GENOME_DIR":"$CONTAINER_GENOME_DIR" \
                "$IMAGE_NAME" \
                bash "$CONTAINER_SCRIPTS_DIR/04_Merging_CTGA.sh" "$@"
            ;;
        
        5|"rider"|"rider-processing")
            echo "Launching 05_Rider_processing.sh script..."
            # USAGE: $0 5 <FLAGS: --merge | --single> [THRESHOLDS/nothing]
            # Examples:
            # nohup ./00_Pipeline_launch.sh 5 --merge 5 > nohup_05_Rider_processing_merge5.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 5 --merge 5 3-5 > nohup_05_Rider_processing_merge5_3-5.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 5 --merge 5 6 7 > nohup_05_Rider_processing_merge5_67.out 2>&1 &

            # Validate arguments
            if [[ "$1" != "--single" && "$1" != "--merge" ]] || [[ -z "$2" || "$2" == --* ]]; then
                echo "[ERROR] Step $STEP requires: <--single|--merge> <min_count> [thresholds]"
                echo "[INFO] You provided: Mode='$1', Min_count='$2'"
                echo "[INFO] Example: $0 $STEP --merge 5 3-5"
                echo -e "\n[INFO] Showing help:\n"
                show_help
                exit 1
            fi

            # Docker execution
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "$PROJECT_DIR":"$CONTAINER_PMAT_DIR" \
                -v "$SCRIPTS_DIR":"$CONTAINER_SCRIPTS_DIR" \
                -v "$GENOME_DIR":"$CONTAINER_GENOME_DIR" \
                -v "$RIDER_DIR":"$CONTAINER_RIDER_DIR" \
                "$IMAGE_NAME" \
                bash "$CONTAINER_SCRIPTS_DIR/05_Rider_processing.sh" "$@"
            ;;

        *)
            echo "[WARNING] Invalid option: $STEP"
            echo "Usage: $0 [1|2|3|4|5] [args...]"
            ;;
    esac
done
