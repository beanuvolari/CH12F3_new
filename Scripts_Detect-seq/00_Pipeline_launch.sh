#!/bin/bash

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

# Save original arguments
ORIGINAL_ARGS=("$@")

for cell in "${CELL_LINE[@]}"; do
  echo "Processing cell line: $cell"

  # Restore original arguments at the start of each iteration
  set -- "${ORIGINAL_ARGS[@]}"

  #exp_dir="${BASE_DIR}/Results/${cell}"
  #mkdir -p "${exp_dir}"

    # Host paths
    #FASTQ_DIR="${BASE_DIR}/Original_fastq/${cell}/merged_fastq"
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


    # Prompt for step if not provided
    if [ -z "$1" ]; then
        echo "Which step do you want to run?"
        echo "1) Fastq to PMAT (01_pmat_generation.sh)"
        echo "2) Filtering (02_Filtering_pmat.sh)"
        echo "3) Sum Mut Loci (03_Sum_mut_loci.sh)"
        echo "4) Merging CTGA (04_Merging_CTGA.sh)"
        echo "5) Rider processing (05_Rider_processing.sh)"
        read -p "Enter the number (0/1/2/3/4): " STEP
    else
        STEP="$1"
    fi

    # Shift the step argument to process the remaining as parameters
    shift

    # Based on selected step
    case "$STEP" in

        1|"fastq-to-pmat")

            # nohup ./00_Pipeline_launch.sh 1 --all > nohup_01_all_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 1 > nohup_01_all_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 1 --single sample_R1.fastq.gz sample2_R1.fastq.gz > nohup_01_single_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 1 --single AID_KO_CIT_merged_R1.fastq.gz > nohup_01_single_AID_KO_CIT_pmat.out 2>&1 &

            # Default behavior: if not specified elaborate all samples in raw_fastq
            MODE="--all"
            SAMPLE_ARGS=()

            # Options parsing for --single or --all
            while (( "$#" )); do
                case "$1" in
                    --single)
                        MODE="--single"
                        shift
                        # Collect files names if they are present
                        if [ $# -eq 0 ]; then
                            # --single option but no sample name provided: interactive request
                            read -p "Provide the name of sample files space-sepaated (es: sample_R1.fastq sample2_R1.fastq): " -a SAMPLE_ARGS
                        else
                            # Save all remaining arguments as sample names
                            SAMPLE_ARGS=("$@")
                        fi
                        break # Exit the while loop after processing --single option
                        ;;
                    --all)
                        MODE="--all"
                        shift
                        ;;
                    *)
                        # Ignore unknown parameters
                        shift
                        ;;
                esac
            done

            echo "Launching 01_pmat_generation.sh script..."

            if [ "$MODE" = "--single" ]; then
                # Single sample mode
                docker run --rm -i \
                    --user $(id -u):$(id -g) \
                    -v "${GENOME_DIR}":"${CONTAINER_GENOME_DIR}" \
                    -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                    -v "${FASTQ_DIR}":"${CONTAINER_FASTQ_DIR}" \
                    -v "${exp_dir}":/scratch \
                    "$IMAGE_NAME2" \
                    bash "$CONTAINER_SCRIPTS_DIR/00_pmat_generation.sh" "$ADAPT1" "$ADAPT2" "--single" "${SAMPLE_ARGS[@]}"
            else
                # Default or --all: analyse all the files in raw_fastq
                docker run --rm -i \
                    --user $(id -u):$(id -g) \
                    -v "${GENOME_DIR}":"${CONTAINER_GENOME_DIR}" \
                    -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                    -v "${FASTQ_DIR}":"${CONTAINER_FASTQ_DIR}" \
                    -v "${exp_dir}":/scratch \
                    "$IMAGE_NAME2" \
                    bash "$CONTAINER_SCRIPTS_DIR/01_pmat_generation.sh" "$ADAPT1" "$ADAPT2" "--all"
            fi
            ;;

        2|"filtering")
            # nohup ./00_Pipeline_launch.sh 2 3-5 > nohup_02_3-5_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 6-10 > nohup_02_6-10_pmat.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 --part2 6 7 8 > nohup_02_678_part2.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 --part2 6-10 > nohup_02_6-10_part2.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 2 --part1 > nohup_02_part1.out 2>&1 &

            echo "Launching 02_Filtering_pmat.sh script..."
            if [ $# -lt 1 ]; then
                read -p "Provide mode flag (e.g., --part1/--part2/--all) and thresholds (space-separated or as range, e.g. 3-6): " -a ARGS
            else
                ARGS=("$@")  
            fi
        
            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "${PROJECT_DIR}":"${CONTAINER_PMAT_DIR}" \
                -v "${SCRIPTS_DIR}":"${CONTAINER_SCRIPTS_DIR}" \
                "${IMAGE_NAME}" \
                bash "${CONTAINER_SCRIPTS_DIR}/02_Filtering_pmat.sh" "${ARGS[@]}"
            ;;

        3|"sum"|"sum-mut-loci")

            # nohup ./00_Pipeline_launch.sh 3 6 > nohup_03_Sum_mut_loci_6.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 3 3 > nohup_03_Sum_mut_loci_3.out 2>&1 &

            echo "Launching 03_Sum_mut_loci.sh script..."
            if [ $# -lt 1 ]; then
                read -p "Provide the threshold from which to start the graphs: " start_threshold
            else
                start_threshold="$1"
            fi

            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "$PROJECT_DIR":"$CONTAINER_PMAT_DIR" \
                -v "$SCRIPTS_DIR":"$CONTAINER_SCRIPTS_DIR" \
                "$IMAGE_NAME" \
                bash "$CONTAINER_SCRIPTS_DIR/03_Sum_mut_loci.sh" "$start_threshold"
            ;;

        4|"merge"|"merging-ctga")

            # nohup ./00_Pipeline_launch.sh 4 > nohup_04_Merging_CTGA.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 4 3-5 > nohup_04_Merging_CTGA_3-5.out 2>&1 &

            echo "Launching 04_Merging_CTGA.sh script..."

            if [ $# -ge 1 ]; then
                THRESHOLDS=("$@")
            else
                THRESHOLDS=()
            fi

            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "$PROJECT_DIR":"$CONTAINER_PMAT_DIR" \
                -v "$SCRIPTS_DIR":"$CONTAINER_SCRIPTS_DIR" \
                -v "$GENOME_DIR":"$CONTAINER_GENOME_DIR" \
                "$IMAGE_NAME" \
                bash "$CONTAINER_SCRIPTS_DIR/04_Merging_CTGA.sh" "${THRESHOLDS[@]}"
            ;;
        
        5|"rider"|"rider-processing")

            # nohup ./00_Pipeline_launch.sh 5 --merge 5 3-5 > nohup_05_Rider_processing_merge5_3-5.out 2>&1 &
            # nohup ./00_Pipeline_launch.sh 5 --merge 5 6 7 > nohup_05_Rider_processing_merge5_67.out 2>&1 &

            echo "Launching 05_Rider_processing.sh script..."

            if [ $# -lt 2 ]; then
                echo "[INFO] You must specify mode (--single or --merge), min_count, and thresholds."
                read -p "Enter mode (--single or --merge): " MODE
                read -p "Enter min_count (e.g. 5): " MIN_COUNT
                read -p "Enter thresholds or range (e.g. 3-6 or 3 5 7): " -a THRESHOLDS
                ARGS=("$MODE" "$MIN_COUNT" "${THRESHOLDS[@]}")
            else
                ARGS=("$@")
            fi

            docker run --rm -i \
                --user $(id -u):$(id -g) \
                -v "$PROJECT_DIR":"$CONTAINER_PMAT_DIR" \
                -v "$SCRIPTS_DIR":"$CONTAINER_SCRIPTS_DIR" \
                -v "$GENOME_DIR":"$CONTAINER_GENOME_DIR" \
                -v "$RIDER_DIR":"$CONTAINER_RIDER_DIR" \
                "$IMAGE_NAME" \
                bash "$CONTAINER_SCRIPTS_DIR/05_Rider_processing.sh" "${ARGS[@]}"
            ;;


        *)
            echo "Invalid option: $STEP"
            echo "Usage: $0 [0|1|2|3|4] [args...]"
            ;;
    esac
done
