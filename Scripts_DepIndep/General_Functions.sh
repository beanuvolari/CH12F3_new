#!/bin/bash

###########################################
#              FUNCTIONS                  #
###########################################

# Functions to record time taken to process files
start_timer() {
    local raw_name="$1"
    local safe_name=$(echo "$raw_name" | sed 's/[^a-zA-Z0-9_]/_/g')

    eval "${safe_name}_start=\$(date +%s)"
    eval "${safe_name}_start_human=\$(date '+%Y-%m-%d %H:%M:%S')"

    local start_human=$(eval echo \${${safe_name}_start_human})
    echo -e "[$raw_name] Started at: $start_human"
}

end_timer() {
    local raw_name="$1"
    local safe_name=$(echo "$raw_name" | sed 's/[^a-zA-Z0-9_]/_/g')

    eval "${safe_name}_end=\$(date +%s)"
    eval "${safe_name}_end_human=\$(date '+%Y-%m-%d %H:%M:%S')"

    local start=$(eval echo \${${safe_name}_start})
    local end=$(eval echo \${${safe_name}_end})
    local end_human=$(eval echo \${${safe_name}_end_human})

    echo -e "[$raw_name] Ended at: $end_human"

    local elapsed=$((end - start))
    local hours=$((elapsed / 3600))
    local minutes=$(( (elapsed % 3600) / 60 ))
    local seconds=$((elapsed % 60))
    echo -e "[$raw_name] Time taken: ${hours}h ${minutes}m ${seconds}s\n"
}


# Function to check if directories already exist or if they need to be created
check_directory_exists() {
    local dir="$1"
    
    if [ -d "$dir" ]; then
        echo -e "[INFO] Directory '$dir' already exists."
    else
        echo -e "\n[WARNING] Directory '$dir' not found!"
        echo -e "[INFO] Creating directory '$dir'..."
        mkdir -p "$dir"
    fi
}

# Function to check if input files exist
check_file_exists() {
    local file="$1"

    if [ -f "$file" ]; then
        echo -e "\n[INFO] Input file exist: $file"
    else
        echo -e "\n[ERROR] File '$file' not found!" >&2
        exit 1
    fi
}

# Function to check if the output file has been created by the script
check_output_generation() {
    local file="$1"
     if [[ ! -f "$file" ]]; then
        echo -e "[ERROR] Output file not found: $file" >&2
        exit 1
    elif [[ ! -s "$file" ]]; then
        echo -e "[ERROR] Output file is empty: $file" >&2
        exit 1
    else
        echo -e "[SUCCESS] Output file exists and is non-empty: $file"
    fi
}

# Function to clear files if they may already exist
clear_file() {
  local file="$1"
  if [ -f "$file" ]; then
    : > "$file"
  fi
}

# Function to remove files that may already exist
remove_file_if_exists() {
    local file="$1"
    if [ -f "$file" ]; then
        rm "$file"
    fi
}

# Function to parses a range like "3-6" or a list like "3 5 7"
parse_thresholds() {
  local input="$1"
  local result=()

  if [[ "$input" =~ ^[0-9]+-[0-9]+$ ]]; then
    IFS="-" read -r start end <<< "$input"
    for ((i=start; i<=end; i++)); do
      result+=("$i")
    done
  else
    result=($input)  # Assume space-separated list
  fi

  echo "${result[@]}"
}

# Simple logging function with timestamp
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}