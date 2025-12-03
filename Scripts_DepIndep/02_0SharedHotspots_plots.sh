#!/bin/bash

# nohup ./02_AID_sharedhotspots_plots.sh CH12F3 3 5 0 -e -u > nohup_sharedhotspots_plots350_eu.out 2>&1 &

# --- CONFIGURATION ---
# Names of the R scripts (must be in the same directory or provide full path)
SCRIPT_DIR="/scratch/scripts"
SCRIPT_EULER="${SCRIPT_DIR}/02_Euler_plot.R"
SCRIPT_VENN="${SCRIPT_DIR}/02_Venn_plot.R"   
SCRIPT_UPSET="${SCRIPT_DIR}/02_Upset_plot.R" 


# --- LOGGING FUNCTIONS ---
# Helper functions to print formatted messages
log_info() { echo -e "[INFO] $1"; }
log_error() { echo -e "[ERROR] $1"; }
log_warn() { echo -e "[WARNING] $1"; }

# --- HELP / USAGE ---
usage() {
    echo "Usage: $0 <cell_line> <threshold> <min_count> <enlargement> [OPTIONS]"
    echo ""
    echo "Positional Arguments (Required in this order):"
    echo "  0. cell_line          Cell line label (e.g., CH12F3)"
    echo "  1. threshold      Numeric value (e.g., 3)"
    echo "  2. min_count      Numeric value (e.g., 5)"
    echo "  3. enlargement    Numeric value (e.g., 1)"
    echo ""
    echo "Options (if none specified, runs ALL by default):"
    echo "  ALL             Runs all scripts (default behavior)"
    echo "  -e, --euler     Runs only Euler.R"
    echo "  -v, --venn      Runs only Venn.R"
    echo "  -u, --upset     Runs only Upset.R"
    echo "  -h, --help      Shows this help message"
    echo ""
    echo "Examples:"
    echo "  $0 CH12F3 3 5 1                    # Runs everything for T=3, MC=5, E=1"
    echo "  $0 CH12F3 3 5 1 ALL                # Explicitly runs everything"
    echo "  $0 CH12F3 3 5 1 -e -u              # Runs only Euler and Upset"
    exit 1
}

# --- ARGUMENT PARSING ---

# State variables
DO_EULER=false
DO_VENN=false
DO_UPSET=false
EXPLICIT_SELECTION=false

# 1. Parse Positional Arguments
# We expect at least 3 arguments (Threshold, MinCount, Enlargement)
if [ $# -lt 4 ]; then
    log_error "Missing required arguments. You must provide cell_line, threshold, min_Count, and enlargement."
    usage
fi

CELL_LINE="$1"
THRESHOLD="$2"
MIN_COUNT="$3"
ENLARGEMENT="$4"

# Shift the first 3 arguments away so we can process flags
shift 4

# Loop through arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -e|--euler)
            DO_EULER=true; EXPLICIT_SELECTION=true; shift ;;
        -v|--venn)
            DO_VENN=true; EXPLICIT_SELECTION=true; shift ;;
        -u|--upset)
            DO_UPSET=true; EXPLICIT_SELECTION=true; shift ;;
        ALL)
            # ALL is implicit if no flags are chosen, but we handle it if typed
            DO_EULER=true
            DO_VENN=true
            DO_UPSET=true
            EXPLICIT_SELECTION=true
            shift # past argument
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_warn "Unknown option ignored: $1"
            shift
            ;;
    esac
done

# --- VALIDATION ---

# Basic validation to ensure inputs are not empty
if [ -z "$CELL_LINE" ] || [ -z "$THRESHOLD" ] || [ -z "$THRESHOLD" ] || [ -z "$ENLARGEMENT" ]; then
    log_error "One of the required parameters (Threshold, MinCount, Enlargement) is empty."
    usage
fi

# Default Logic: If the user didn't specify any flags, enable ALL
if [ "$EXPLICIT_SELECTION" = false ]; then
    log_info "No specific flags detected. Running mode: ALL"
    DO_EULER=true
    DO_VENN=true
    DO_UPSET=true
fi

# --- EXECUTION ---

log_info "Starting pipeline for $CELL_LINE -> T=${THRESHOLD}, MC=${MIN_COUNT}, E=${ENLARGEMENT}"

# Define input file path based on parameters
BASE_DIR="/scratch/DepIndep_dataset/Threshold_${THRESHOLD}/Enlargement_${ENLARGEMENT}"
input_file="${BASE_DIR}/${ENLARGEMENT}_${CELL_LINE}_SharedHotspots_${THRESHOLD}_rider${MIN_COUNT}.bed"
OUTPUT_DIR="${BASE_DIR}/Plots_Euler_Venn_Upset"
[ -d "$OUTPUT_DIR" ] || { echo "[INFO] Output_dir doesn't exist. Creating it..." ; mkdir -p "$OUTPUT_DIR"; }


# --- INPUT FILE VALIDATION ---
log_info "Checking that the dataset is fully tab-delimited..."

# Controlla se esiste almeno una riga senza tabulazioni
if awk 'index($0, "\t") == 0 { print NR; exit 1 }' "$input_file"; then
    # Se awk termina con exit 0, significa che NON ha trovato linee errate
    log_info "Input file is correctly tab-delimited."
else
    # awk ha trovato almeno una riga senza tab → errore
    echo "[ERROR] Your input file shows formatting problems."
    echo "[ERROR] The script needs a tab-separated dataset."
    echo "[ERROR] First problematic line number shown above."
    exit 1
fi

# If input file has problems, uncomment the following lines to fix common issues:
# 1. Remove Windows carriage returns (\r) if any
#sed -i 's/\r//g' "$input_file"

# 2. Remove ALL trailing tabs and spaces at the end of lines
# This regex matches any combination of tab/space at the end ($) and deletes it
#sed -i 's/[[:space:]]*$//' "$input_file"

# Define a function to run R scripts to avoid code repetition
run_r_script() {
    local script_path="$1"
    local input_file="$2"
    local output_dir="$3"
    local script_name=$(basename "$script_path")

    if [ -f "$script_path" ]; then
        log_info "Launching $script_name..."
        # We pass ALL 3 arguments to R: Threshold, MinCount, Enlargement
        Rscript "$script_path" "$input_file" "$output_dir" "$THRESHOLD" "$MIN_COUNT" "$ENLARGEMENT"
        
        if [ $? -eq 0 ]; then 
            log_info "$script_name completed successfully."
        else 
            log_error "$script_name failed."
        fi
    else
        log_error "Script file not found: $script_path"
    fi
}

# 1. Execute EULER Script
if [ "$DO_EULER" = true ]; then
    run_r_script "$SCRIPT_EULER" "$input_file" "$OUTPUT_DIR"
fi

# 2. Execute VENN Script
if [ "$DO_VENN" = true ]; then
    run_r_script "$SCRIPT_VENN" "$input_file" "$OUTPUT_DIR"
fi

# 3. Execute UPSET Script
if [ "$DO_UPSET" = true ]; then
    run_r_script "$SCRIPT_UPSET" "$input_file" "$OUTPUT_DIR"
fi

log_info "All requested operations finished."