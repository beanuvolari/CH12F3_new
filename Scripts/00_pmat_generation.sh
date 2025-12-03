#!/bin/bash


# Variables
raw_dir="/scratch/raw_fastq"
genome_dir="/scratch/genome"

# Output directories creation
fix_fastq_dir="/scratch/fix_fastq"
bam_dir="/scratch/bam_hisat3n"
pmat_dir="/scratch/pmat_and_mpmat"


##########################
#       FUNCTIONS        #
##########################

# Source general helper functions
source /scratch/scripts/General_Functions.sh

# Check if file exists function
file_exists() {
    if [ -f "$1" ]; then
        return 0
    else
        return 1
    fi
}

# Build BWA index if needed
build_bwa_index() {
    local genome_file=$1 
    local filename_g=$2
    local bwa_index_check="${genome_dir}/bwa_${filename_g}/${filename_g}.fa.bwt"

    if file_exists "${bwa_index_check}"; then
        echo "[INFO] BWA index found. Skipping index generation."
    else
        echo "[INFO] BWA index not found. Generating index..."
        samtools faidx "${genome_file}"

        mkdir -p "${genome_dir}/bwa_${filename_g}"
        cp "${genome_file}"* "${genome_dir}/bwa_${filename_g}/"
        cd "${genome_dir}/bwa_${filename_g}" || { echo "[ERROR] Failed to enter BWA index directory"; exit 1; }
        bwa index "${genome_file}"
        cd - > /dev/null
        echo "[INFO] Current directory after bwa index generation is: $(pwd)"
    fi
}

# Build HISAT3N index if needed
build_hisat3n_index() {
    local genome_file=$1
    local filename_g=$2
    local hisat_index_check="${genome_dir}/hisat3n_${filename_g}_CT/${filename_g}.fa.3n.CT.1.ht2"

    if file_exists "${hisat_index_check}"; then
        echo "[INFO] HISAT3N index found. Skipping index generation."
    else
        echo "[INFO] HISAT3N index not found. Generating index..."
        mkdir -p "${genome_dir}/hisat3n_${filename_g}_CT"
        cp "${genome_dir}/${genome_file}"* "${genome_dir}/hisat3n_${filename_g}_CT"
        cd "${genome_dir}/hisat3n_${filename_g}_CT" || { echo "[ERROR] Failed to enter HISAT3N index directory"; exit 1; }
        /home/hisat-3n/hisat-3n-build --base-change C,T "${genome_file}" "${genome_file}" > "hisat3n_${filename_g}_CT_index.log"
        cd - > /dev/null
        echo "[INFO] Current directory after hisat3n index generation is: $(pwd)"
    fi
}

# Detect genome fasta and build indexes function
detection_fasta() {
    echo "[INFO] Starting genome fasta detection..."

    cd "${genome_dir}" || { echo "[ERROR] Genome directory not found!"; exit 1; }
    
    # Auto-detect genome FASTA file
    genome_file=$(ls *.fa* | head -n 1)
    if [ -z "$genome_file" ]; then
        echo "[ERROR] Genome FASTA file not found in ${genome_dir}!"
        exit 1
    fi
    file_g=$(basename -- "$genome_file")
    filename_g="${file_g%.*}"

    build_bwa_index "$genome_file" "$filename_g"
    build_hisat3n_index "$genome_file" "$filename_g"
}


# Extract file extension and base sample name from FASTQ filename
extract_sample_info() {
    local filename=$1
    local ext=""
    local sample_base=""

    # Detect reliable file extension
    if [[ $filename == *.fastq.gz ]]; then
        ext="fastq.gz"
    elif [[ $filename == *.fq.gz ]]; then
        ext="fq.gz"
    else
        # If not a known extension, assign empty or handle differently if needed
        ext=""
    fi

    # Extract base sample name by removing extension and suffixes _R1 or _1 if present
    local base_with_ext=$(basename "$filename" ".$ext")
    base_with_ext=${base_with_ext%_R1}
    base_with_ext=${base_with_ext%_1}
    sample_base="$base_with_ext"

    # Return extension and base sample name as space-separated strings
    echo "$ext $sample_base"
}

# Extract extension and sample info from single (provided) files
extract_sample_info_single() {
  local files=("$@")
  for fqfile in "${files[@]}"; do
    read ext sample <<< $(extract_sample_info "$fqfile")
    echo "$ext $sample"
  done
}

# Extract extension and sample info from all files in raw_fastq folder
extract_sample_info_all() {
  for fqfile in "${raw_dir}"/*_R1.fastq.gz "${raw_dir}"/*_R1.fq.gz "${raw_dir}"/*_1.fastq.gz "${raw_dir}"/*_1.fq.gz; do
    if [[ -f $fqfile ]]; then
      read ext sample <<< $(extract_sample_info "$fqfile")
      echo "$ext $sample"
    fi
  done
}




#------Processing steps functions------#

# Adapter trimming function
trim_adapters() {
    local sample=$1
    local adapt1=$2
    local adapt2=$3
    local ext=$4

    echo "[INFO] Starting adapter trimming for sample ${sample}..."

    mkdir -p "${fix_fastq_dir}"
    local in_fq_R1="${raw_dir}/${sample}_R1.${ext}"
    local in_fq_R2="${raw_dir}/${sample}_R2.${ext}"
    local out_fq_R1="${fix_fastq_dir}/${sample}_R1_cutadapt.fastq.gz"
    local out_fq_R2="${fix_fastq_dir}/${sample}_R2_cutadapt.fastq.gz"
    local log="${fix_fastq_dir}/${sample}_cutadapt.log"

    cutadapt -j 0 --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 55 -a "${adapt1}" -A "${adapt2}" -o "${out_fq_R1}" -p "${out_fq_R2}" "${in_fq_R1}" "${in_fq_R2}" > "${log}"

    echo "[INFO] Adapter trimming completed for sample ${sample}."
}

# HISAT3N alignment function
align_hisat3n() {
    local sample=$1
    local genome_file=$2
    local filename_g=$3

    echo "[INFO] Starting HISAT3N alignment for sample ${sample}..."

    mkdir -p "${bam_dir}"

    local in_fq_R1="${fix_fastq_dir}/${sample}_R1_cutadapt.fastq.gz"
    local in_fq_R2="${fix_fastq_dir}/${sample}_R2_cutadapt.fastq.gz"
    local out_bam="${bam_dir}/${sample}_hisat3n.bam"
    local unmapped_fq="${bam_dir}/${sample}_hisat3n_unmapped.fastq.gz"
    local log="${bam_dir}/${sample}_hisat3n.log"
    local ref_idx="${genome_dir}/hisat3n_${filename_g}_CT/${genome_file}"

    /home/hisat-3n/hisat-3n -x "${ref_idx}" -1 "${in_fq_R1}" -2 "${in_fq_R2}" -p 20 --sensitive --base-change C,T --unique-only --repeat-limit 1000 --no-spliced-alignment -X 700 --un-conc-gz "${unmapped_fq}" --summary-file "${log}" --rg-id "${sample}" --rg "PL:ILLUMINA" --rg "ID:${sample}" --rg "SM:${sample}" | samtools view -hb > "${out_bam}"

    echo "[INFO] HISAT3N alignment completed for sample ${sample}."
}

# Select reads with low mapping quality
select_low_quality_reads() {
    local sample=$1

    echo "[INFO] Selecting low MAPQ reads for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_hisat3n.bam"
    local out_bam="${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.bam"

    samtools view -h -@ 4 "${in_bam}" | awk '$1 ~ "@" || $5 <= 20 {print $0}' | samtools view -@ 4 -hb > "${out_bam}"

    echo "[INFO] Low MAPQ reads selected for sample ${sample}."
}

# Sort BAM by read name
sort_bam_by_name() {
    local sample=$1

    echo "[INFO] Sorting BAM file by read name for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.bam"
    local out_bam="${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam"
    local temp_file="${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam.temp"

    samtools sort -O BAM -o "${out_bam}" -T "${temp_file}" -@ 15 -m 2G -n "${in_bam}"

    echo "[INFO] BAM sorting by read name completed for sample ${sample}."
}

# Extract low quality reads to FASTQ
extract_low_quality_reads() {
    local sample=$1
    local genome_file=$2
    local filename_g=$3

    echo "[INFO] Extracting low quality reads to FASTQ for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam"
    local ref_genome_fa="${genome_dir}/hisat3n_${filename_g}_CT/${genome_file}"
    local out_fq_R1="${bam_dir}/${sample}_hisat3n.LowerMAPQ20_R1.fastq.gz"
    local out_fq_R2="${bam_dir}/${sample}_hisat3n.LowerMAPQ20_R2.fastq.gz"

    samtools fastq -@ 15 -0 /dev/null -s /dev/null -n -F 0x900 -1 "${out_fq_R1}" -2 "${out_fq_R2}" --reference "${ref_genome_fa}" "${in_bam}"

    echo "[INFO] Extraction to FASTQ completed for sample ${sample}."
}

# Merge unmapped and low quality FASTQ reads
merge_unmapped_and_low_quality_reads() {
    local sample=$1

    echo "[INFO] Merging unmapped and low quality reads for sample ${sample}..."

    local low_fq_R1="${bam_dir}/${sample}_hisat3n.LowerMAPQ20_R1.fastq.gz"
    local low_fq_R2="${bam_dir}/${sample}_hisat3n.LowerMAPQ20_R2.fastq.gz"
    local unmapped_fq_R1="${bam_dir}/${sample}_hisat3n_unmapped.fastq.1.gz"
    local unmapped_fq_R2="${bam_dir}/${sample}_hisat3n_unmapped.fastq.2.gz"
    local out_fq_R1="${bam_dir}/${sample}_R1_unmapped_and_LowerMAPQ20.fastq.gz"
    local out_fq_R2="${bam_dir}/${sample}_R2_unmapped_and_LowerMAPQ20.fastq.gz"

    cat "${low_fq_R1}" "${unmapped_fq_R1}" > "${out_fq_R1}"
    cat "${low_fq_R2}" "${unmapped_fq_R2}" > "${out_fq_R2}"

    echo "[INFO] Merging completed for sample ${sample}."
}

# Realign with BWA MEM
realign_bwa_mem() {
    local sample=$1
    local genome_file=$2
    local filename_g=$3

    echo "[INFO] Starting BWA MEM realignment for sample ${sample}..."

    local in_fq_R1="${bam_dir}/${sample}_R1_unmapped_and_LowerMAPQ20.fastq.gz"
    local in_fq_R2="${bam_dir}/${sample}_R2_unmapped_and_LowerMAPQ20.fastq.gz"
    local bwa_index="${genome_dir}/bwa_${filename_g}/${genome_file}"
    local out_bam="${bam_dir}/${sample}_bwa_realign.bam"
    local bwa_log="${bam_dir}/${sample}_bwa_realign.log"

    bwa mem "${bwa_index}" "${in_fq_R1}" "${in_fq_R2}" -t 20 -M -R "@RG\tID:${sample}\tPL:ILLUMINA\tSM:${sample}" 2>"${bwa_log}" | samtools view -h -b -q 20 -f 3 -F 256 > "${out_bam}"

    echo "[INFO] BWA MEM realignment completed for sample ${sample}."
}

# Merge HISAT3N and BWA BAM files
merge_bam_files() {
    local sample=$1

    echo "[INFO] Merging HISAT3N and BWA BAM files for sample ${sample}..."

    local in_bam_bwa="${bam_dir}/${sample}_bwa_realign.bam"
    local in_bam_hisat3n="${bam_dir}/${sample}_hisat3n.bam"
    local out_bam="${bam_dir}/${sample}_merge.MAPQ20.bam"

    samtools cat -o "${out_bam}" "${in_bam_hisat3n}" "${in_bam_bwa}"

    echo "[INFO] BAM files merged for sample ${sample}."
}

# Final BAM sort by genomic coordinates
final_sort_bam() {
    local sample=$1

    echo "[INFO] Sorting merged BAM file by genomic coordinates for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_merge.MAPQ20.bam"
    local out_bam="${bam_dir}/${sample}_merge_sort.MAPQ20.bam"
    local temp_file="${bam_dir}/${sample}_merge_sort.MAPQ20.bam.temp"

    samtools sort -O BAM -o "${out_bam}" -T "${temp_file}" -@ 15 -m 2G "${in_bam}"

    echo "[INFO] Sorting completed for sample ${sample}."
}

# Remove duplicates with Picard
remove_duplicates() {
    local sample=$1

    echo "[INFO] Removing duplicates for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_merge_sort.MAPQ20.bam"
    local out_bam="${bam_dir}/${sample}_merge_sort_rmdup.MAPQ20.WithClip.bam"
    local out_log="${bam_dir}/${sample}_merge_sort_rmdup.log"
    local out_matrix="${bam_dir}/${sample}_merge_sort_rmdup.matrix"

    java -Xms50g -Xmx50g -XX:ParallelGCThreads=20 -jar /home/picard.jar MarkDuplicates I="${in_bam}" O="${out_bam}" M="${out_matrix}" ASO=coordinate REMOVE_DUPLICATES=true 2>"${out_log}"

    echo "[INFO] Duplicate removal completed for sample ${sample}."
}

# Filter BAM for clipping, concordancy, MAPQ, and secondary alignments
filter_bam() {
    local sample=$1
    local genome_file=$2
    local filename_g=$3

    echo "[INFO] Filtering BAM file for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_merge_sort_rmdup.MAPQ20.WithClip.bam"
    local out_bam="${bam_dir}/${sample}_merge_sort_rmdup.MAPQ20.bam"
    local ref_genome_fa="${genome_dir}/hisat3n_${filename_g}_CT/${genome_file}"

    samtools view -@ 4 -h "${in_bam}" -q 20 -f 3 -F 256 | /home/samclip --ref "${ref_genome_fa}" --max 3 --progress 0 | \
    awk 'function abs(v) {return v < 0 ? -v : v} $1 ~ "@" || ($7 == "=" && abs($9) <= 2500) {print $0}' | samtools view -@ 4 -hb > "${out_bam}"

    echo "[INFO] BAM filtering completed for sample ${sample}."
}

# Create BAM index
index_bam() {
    local sample=$1

    echo "[INFO] Creating BAM index for sample ${sample}..."

    local in_bam="${bam_dir}/${sample}_merge_sort_rmdup.MAPQ20.bam"
    local out_bam_index="${in_bam}.bai"

    samtools index -@ 10 "${in_bam}" "${out_bam_index}"

    echo "[INFO] BAM index created for sample ${sample}."
}

# Convert BAM to PMAT format
convert_to_pmat() {
    local sample=$1
    local genome_file=$2
    local filename_g=$3

    echo "[INFO] Converting BAM to PMAT for sample ${sample}..."

    mkdir -p "${pmat_dir}"

    local in_bam="${bam_dir}/${sample}_merge_sort_rmdup.MAPQ20.bam"
    local ref_genome_fa="${genome_dir}/hisat3n_${filename_g}_CT/${genome_file}"
    local out_pmat="${pmat_dir}/${sample}_CLEAN.pmat"
    local out_log="${pmat_dir}/${sample}_CLEAN.log"

    /root/anaconda3/envs/DetectSeq/bin/python /home/Detect-seq/src/detect_seq/bam2pmat.py -i "${in_bam}" -r "${ref_genome_fa}" -o "${out_pmat}" -p 20 --out_format pmat --bed_like_format True --mut_type ALL --block_size 100000 --cover_num_cutoff 0 --mut_num_cutoff 0 --mut_ratio_cutoff 0 --keep_temp_file False --out_header False > "${out_log}"

    echo "[INFO] PMAT conversion completed for sample ${sample}."
}

# Cleanup temporary files
cleanup_temp_files() {
    local sample=$1

    echo "[INFO] Cleaning up temporary files for sample ${sample}..."

    rm -f "${fix_fastq_dir}/${sample}_R1_cutadapt.fastq.gz"
    rm -f "${fix_fastq_dir}/${sample}_R2_cutadapt.fastq.gz"
    rm -f "${bam_dir}/${sample}_hisat3n.bam"
    rm -f "${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.bam"
    rm -f "${bam_dir}/${sample}_hisat3n.MAPQ20.LowerMAPQ20.SortName.bam"
    rm -f "${bam_dir}/${sample}_hisat3n.LowerMAPQ20_R1.fastq.gz"
    rm -f "${bam_dir}/${sample}_hisat3n.LowerMAPQ20_R2.fastq.gz"
    rm -f "${bam_dir}/${sample}_hisat3n_unmapped.fastq.1.gz"
    rm -f "${bam_dir}/${sample}_hisat3n_unmapped.fastq.2.gz"
    rm -f "${bam_dir}/${sample}_R1_unmapped_and_LowerMAPQ20.fastq.gz"
    rm -f "${bam_dir}/${sample}_R2_unmapped_and_LowerMAPQ20.fastq.gz"
    rm -f "${bam_dir}/${sample}_bwa_realign.bam"
    rm -f "${bam_dir}/${sample}_merge.MAPQ20.bam"
    rm -f "${bam_dir}/${sample}_merge_sort.MAPQ20.bam"
    rm -f "${bam_dir}/${sample}_merge_sort_rmdup.MAPQ20.WithClip.bam"

    echo "[INFO] Cleanup completed for sample ${sample}."
}


#------Main processing functions-------#

process_sample() {
    local sample=$1
    local adapt1=$2
    local adapt2=$3
    local genome_file=$4
    local filename_g=$5
    local ext=$6

    # Start timing
    start_timer "Analysis_${sample}"
    trim_adapters "$sample" "$adapt1" "$adapt2" "$ext"
    align_hisat3n "$sample" "$genome_file" "$filename_g"
    select_low_quality_reads "$sample"
    sort_bam_by_name "$sample"
    extract_low_quality_reads "$sample" "$genome_file" "$filename_g"
    merge_unmapped_and_low_quality_reads "$sample"
    realign_bwa_mem "$sample" "$genome_file" "$filename_g"
    merge_bam_files "$sample"
    final_sort_bam "$sample"
    remove_duplicates "$sample"
    filter_bam "$sample" "$genome_file" "$filename_g"
    index_bam "$sample"
    convert_to_pmat "$sample" "$genome_file" "$filename_g"
    cleanup_temp_files "$sample"
    echo -e "[INFO] Processing of sample ${sample} completed.\n"
    end_timer "Analysis_${sample}"
    # End timing
}

# Function to process single samples specified by file names
detect_seq_single() {
  local adapt1=$1
  local adapt2=$2
  shift 2
  local sample_files=("$@")

  extract_sample_info_single "${sample_files[@]}" | while read ext sample; do
    process_sample "$sample" "$adapt1" "$adapt2" "$genome_file" "$filename_g" "$ext"
  done
}

# Function to process all samples found in raw_fastq dir
detect_seq_all() {
  local adapt1=$1
  local adapt2=$2

  extract_sample_info_all | while read ext sample; do
    process_sample "$sample" "$adapt1" "$adapt2" "$genome_file" "$filename_g" "$ext"
  done
}


##########################
#          MAIN          #
##########################

# Parsing arguments and setting mode
start_timer "Total_Analysis"

# Check minimum presence of adapters
if [ "$#" -lt 2 ]; then
    echo "[ERROR] Usage: $0 adapt1 adapt2 [--single sample1 sample2 ... | --all]"
    exit 1
fi

adapt1="$1"
adapt2="$2"
shift 2

MODE="all"
samples_to_process=""

if [ "$#" -gt 0 ]; then
    case "$1" in
        --single)
            MODE="single"
            shift
            if [ "$#" -eq 0 ]; then
                # No samples provided after --single, request interactively
                 echo "[WARNING] No samples specified for --single mode."
                read -p "Enter sample file names space-separated: " -a interactive_samples
                samples_to_process=$(printf '%s,' "${interactive_samples[@]}")
                samples_to_process=${samples_to_process%,}
            else
                # Take all samples provided after --single
                samples_to_process=$(printf '%s,' "$@")
                samples_to_process=${samples_to_process%,}
            fi
            ;;
        --all)
            MODE="all"
            shift
            ;;
        *)
            echo "[WARNING] Unknown parameter: $1, ignored."
            ;;
    esac
fi


# Run genome fasta detection and indexing
detection_fasta

# Choose which processing function to run based on mode
if [ "$MODE" = "single" ]; then
  echo "[INFO] Running in single mode for samples: $samples_to_process"
  IFS=',' read -r -a sample_array <<< "$samples_to_process"
  detect_seq_single "$adapt1" "$adapt2" "${sample_array[@]}"
else
  echo "[INFO] Running in all mode; processing all samples in raw_fastq directory"
  detect_seq_all "$adapt1" "$adapt2"
fi

end_timer "Total_Analysis"
echo "[INFO] All processing completed."
