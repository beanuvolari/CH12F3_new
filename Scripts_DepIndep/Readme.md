# Define AID dependent and independent hotspots

## 00_DepIndep_launch

This repository contains a pipeline to define and visualize AID-dependent and independent hotspots. The main script `00_DepIndep_launch.sh` acts as a launcher to orchestrate the different steps of the analysis using Docker containers. This ensures modularity and reproducibility across different environments.

***Prerequisites***

- Docker must be installed and running.
- The Docker image detectseqpipe:4 (or the version specified in the script).
```
docker pull detectseqpipe:4
```

***Configuration***

To run this script, update the following variables at the top of the script: 

- BASE_DIR: the **absolute paths** to the project folder;
```bash
# Base directory for the project
BASE_DIR="Path-to-your-folder/Detect-seq_Project/"
```
- ARGUMENTS:
    - CELL_LINE: **name of the cell line** to process;
    - THRESHOLD: minimum n° of mutations to consider **significant** a **mutation position**
    - MIN_COUNT: minimum n° of mutations to consider **significant** a **hotspot**
    - ENLARGEMENT: Genomic **window expansion**.  
        - 0 = No enlargement
        - 1 = ±1000 bp  (start-1000bp; end+1000bp)
        - 2 = ±2000 bp  (start-2000bp; end+2000bp)
        - 3 = ±4000 bp  (start-4000bp; end+4000bp)
```bash
# Cell lines to process
CELL_LINE="cell_line1"
THRESHOLD="threshold"
MIN_COUNT="min_count"
ENLARGEMENT="enlargement"
```
- ORDERED_SAMPLES: array with the sample names to process, in the order in which they should be analyzed.
```bash
# Samples'order
ORDERED_SAMPLES=("sample_1" "sample_2" "sample_3" ...)
```


***Usage***  
```bash
# Syntax
./00_DepIndep_launch.sh <step> [options...]

# Help
./00_DepIndep_launch.sh --help
./00_DepIndep_launch.sh -h
```
- If no 'step' is specified, the script will prompt you an error. Have a look at the help.  
- The script automatically handles volume mapping between your host directories and the Docker container (/scratch/...)


## Steps


### **1) Generation of a cross-sample shared hotspots dataset**

Option: `1` or `Shared_Hotspots`

This step runs `01_0SharedHotspots.sh`. It wraps `01_1SharedHotspots.R` that collects and merges hotspot regions ("RIDER" regions) detected across multiple samples, producing a consolidated dataset. It allows for the optional expansion of hotspot coordinates via the **enlargement** parameter.


***Input Requirement:*** The script expects the RIDER output:
- Location: `/scratch/pmat/Filtered_pmat/<SAMPLE>/merged_CTGA/threshold_<T>/<T>_rider<MC>`  
- Filename format: `<sample_name>_merged_CTGA_<T>_sorted-RIDER.clean*.bed`


***Usage:***
```bash
# Run:
./00_DepIndep_launch.sh 1
```
Batch mode (recommended)
```bash
# Run via command line (parameters are defined above in the script config part)

# Syntax:
nohup ./00_DepIndep_launch.sh 1 > nohup_01SharedHotspots_merge.out 2>&1 &
```
Example:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 1 > nohup_01SharedHotspots_merge350.out 2>&1 &
```


***Output:***  
The output table is saved in 
- Location: `DepIndep_dataset/Threshold_<T>/Enlargement_<E>/`

- Filename format: `<E>_<CELL>_SharedHotspots_<T>_rider<MC>.bed` 

The script is designed to be run in a **batch mode**.  

___

### **2) Plotting (Euler, Upset, Venn)**

Option: `2` or `Plots`

This step runs `02_0SharedHotspots_plots.sh` (which wraps three R scripts: `02_Euler_plot.R`, `02_Upset_plot.R` and `02_Venn_plot.R`). It generates visualizations to analyze the intersections of the hotspots defined in Step 1.  

***Input Requirement:*** The script expects a tab-delimited BED-like file generated in Step 1.

- Location: `DepIndep_dataset/Threshold_<T>/Enlargement_<E>`

- Filename format: `<E>_<CELL>_SharedHotspots_<T>_rider<MC>.bed` 


***Usage & Flags:***

You can run all plots at once or select specific ones using flags.

```bash
# Syntax:
./00_DepIndep_launch.sh 2 [FLAGS: -v (venn), -e (euler), -u (upset)]
```

```
Flag            Description
ALL             Default. Runs Euler, Venn, and UpSet plots.
-e or --euler   Runs only the Euler plot (proportional overlaps).
-v or --venn    Runs only the Venn diagram (limit: max 5 sets).
-u or --upset   Runs only the UpSet plot (complex intersections).
```
  
Batch mode (recommended)
```bash
# Run via command line (parameters are defined above in the script config prt)

# Generate ALL plots (Default):
nohup ./00_DepIndep_launch.sh > nohup_DepIndep_plot.out 2>&1 &

# Generate only Euler and UpSet plots:
nohup ./00_DepIndep_launch.sh 2 -e -u > nohup_DepIndep_ploteu.out 2>&1 &
```

Examples:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0

nohup ./00_DepIndep_launch.sh > nohup_DepIndep_plot350.out 2>&1 &

nohup ./00_DepIndep_launch.sh 2 -e -u > nohup_DepIndep_ploteu350.out 2>&1 &
```


If you run `./00_DepIndep_launch.sh 2` without flags, the script will generate all the plot types.

***Outputs:*** Plots are saved as high-resolution PNGs (300 DPI) in: `DepIndep_dataset/Threshold_<T>/Enlargement_<E>/Plots_Euler_Venn_Upset/`
- Euler: `<E>_Plots_Euler_<T>_<MC>.png`
- Venn:  `<E>_Plots_Venn_<T>_<MC>.png`
- Upset: `<E>_Plots_Upset_<T>_<MC>.png`


***Note on R Scripts:***

- Venn: will automatically skip generation if the input dataset contains more than 5 sample columns.

- Euler: uses the eulerr package. It may fail if the set intersections are too complex to map geometrically.

- UpSet: converts "+" values to binary (0/1) for intersection analysis.  

***Troubleshooting***  
- `"[ERROR] Your input file shows formatting problems"`: The plotting script checks if the input file is strictly tab-delimited. If this error occurs, ensure your dataset does not contain spaces instead of tabs or malformed line endings (e.g., Windows \r).

- Venn Plot missing: Check the logs. If your data has >5 sample columns, the Venn diagram is intentionally skipped.
___

### **3) AID-Dependent/ AID-Independent hotspots split**

Option: `3` or `AID_DepIndep_Split` 

This step runs the script `03_AID_DepIndephotspots_split.sh`, which separates the shared hotspot dataset (generated in Step 1) into two categories for a user-selected sample:  
- **AID-Dependent** hotspots: regions where the selected sample carries a “**–**”

- **AID-Independent** hotspots: regions where the selected sample carries a “**+**”

The script identifies the correct sample column in the input file and automatically generates two BED-like outputs.

**Input Requirement:**
This step expects as input the shared hotspots file created in Step 1:

Location:
`DepIndep_dataset/Threshold_<T>/Enlargement_<E>/`

Filename format:
`<E>_<CELL>_SharedHotspots_<T>_rider<MC>.bed`


***Parameters:***

You must supply one parameter: 
- REF_SAMPLE_NAME – The specific sample whose “+/–” column will define Dependent vs Independent hotspots.


***Usage:***
```bash
# Syntax:
./00_DepIndep_launch.sh 3 <REF_SAMPLE_NAME>
```

Batch mode (recommended)
```bash
# Run via command line (parameters are defined above in the script config part)

# Syntax:
nohup ./00_DepIndep_launch.sh 3 <REF_SAMPLE_NAME> > nohup_03AID_DepIndephotspots_split.out 2>&1 &
```

Examples:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh AID_KO > nohup_03AID_DepIndephotspots_split350.out 2>&1 &
```

***Outputs:*** Two table are generate
- Dependent:
    - Location:`/scratch/DepIndep_dataset/Threshold_<T>/Enlargement_<E>/Dependent_hotspots/`
    - FIlename: `<E>_<CL>_AID_Dependent_hotspots_<T>_<MC>.bed`
- Independent:
    - Location:`/scratch/DepIndep_dataset/Threshold_<T>/Enlargement_<E>/Independent_hotspots/`
    - FIlename: `<E>_<CL>_AID_Independent_hotspots_<T>_<MC>.bed`
___


### **4) Ranking of AID-Dependent/Independent Hotspots**

Option: `4` or `Ranking_AID_DepIndep`

This step runs `04_0Ranking_AID_DepIndep_hotspots.sh`, which wraps `04_1Ranking_AID_DepIndep_hotspots.R`. It annotates the Dependent or Independent hotspots with the actual mutation counts and density from a specific sample's PMAT data. The script produces two separate rankings for each sample: one based on mutation density and one based on the absolute number of mutations.

***Input Requirement:*** 
- **Hotspots:** The BED files generated in Step 3 (Dependent/Independent).
    - `<E>_<CL>_AID_Dependent_hotspots_<T>_<MC>.bed`
    - `<E>_<CL>_AID_Independent_hotspots_<T>_<MC>.bed`
- **Mutations:** The sorted BED mutation file located in: `/scratch/pmat/Filtered_pmat/<SAMPLE>/merged_CTGA/threshold_<T>/`.

***Usage & Flags:***  
You must provide the name of the sample on which you want base the ranking. Optional flags allow you to restrict the analysis to specific types or ranking modes.

```bash
# Syntax:
./00_DepIndep_launch.sh 4 <ANALYSIS_SAMPLE_NAME> [FLAGS: --Dep / --Indep] [--mut / --dens]
```
```
Flag            Description

--Dep           Process only Dependent file
--Indep         Process only Independent file
none/ALL        Process both Depenedent and Independent files

--mut           Ranks hotspots by number of mutations
--dens          Ranks hotspots by density
none/ALL        Ranks hotspots by both mutations and density
```
Batch mode (recommended)
 
```bash
# Full analysis for a specific sample:
nohup ./00_DepIndep_launch.sh 4 UM_TKO_CIT_Duv > nohup_04Ranking.out 2>&1 &

# Specific analysis (Dependent hotspots ranked by density):
nohup ./00_DepIndep_launch.sh 4 UM_TKO_CIT_Duv --Dep --dens > nohup_04RankingDepdens.out 2>&1 &
```
Examples:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 4 UM_TKO_CIT_Duv --Dep --dens > nohup_04RankingDepdens350.out 2>&1 &
```
***Outputs:***  
 Ranked TSV tables are saved in: 
 - Location: `DepIndep_dataset/Threshold_<T>/Enlargement_<E>/<TYPE>endent_hotspots/<SAMPLE>/Ranking_<TYPE>_hotspots/<MODE>/Full_hotspots/`
 - Filename:   
 `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_<T><MC>.tsv`

 ___

### **5) Filtering Fragile Sites (Blacklist Removal)**
Option: `5` or `Filtering_FragileSites`

This step runs `05_Remove_fragile_sites.R`. It intersects the ranked hotspots with a known blacklist (e.g., mm9-blacklist.bed) to remove regions associated with genomic fragility or sequencing artifacts.

***Input Requirement:***
 - Ranked TSV files from Step 4.
    `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_<T><MC>.tsv`
 - Blacklist file: `/scratch/reference_data/mm9-blacklist.bed`.

***Usage:***
```bash
# Syntax:
./00_DepIndep_launch.sh 5 <ANALYSIS_SAMPLE_NAME> [FLAGS: --Dep / --Indep] [--mut / --dens]
```
Batch mode (recommended)
```bash
# Full analysis for a specific sample:
nohup ./00_DepIndep_launch.sh 5 UM_TKO_CIT_Duv > nohup_05Filtering_FragileSites.out 2>&1 &
```
Example:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 5 UM_TKO_CIT_Duv --Dep --mut > nohup_05Filtering_FragileSites_Depmut350.out 2>&1 &
```

***Outputs:***  
 Files are created in a Filtered_hotspots subdirectory:
 - Filtered: `<...Ranked_hotspots_..._blackfilt.tsv` (hotspots passing the filter).
 - Deleted: `<...Deleted_hotspots_..._blacklist.tsv` (hotspots removed).
 ___


### **6) On/Off-Target Finding**
Option: `6` or `OnOff_DepIndep_finding`

This step runs `06_OnOff_DepIndep_finding.R`. It identifies which hotspots overlap with known HTGTS ON-target and OFF-target regions. The script automatically processes both the "Full" ranked hotspots and the "Filtered" (blackfilt) hotspots generated in the previous steps.

***Input Requirement:***
- Ranked/Filtered TSV files from Steps 4 and 5:
    - Ranked: `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_<T><MC>.tsv`
    - Filtered: `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_<T><MC>_blackfilt.tsv` 
- Reference files: 
    - Location: `/scratch/reference_data/`
        - `ON_target.bed`
        - `OFF_target.bed`


***Usage:***
```bash
# Syntax:
./00_DepIndep_launch.sh 6 <ANALYSIS_SAMPLE_NAME> [FLAGS: --Dep / --Indep] [--mut / --dens]
```
Batch mode (recommended)
```bash
# Full analysis for a specific sample:
nohup ./00_DepIndep_launch.sh 6 UM_TKO_CIT_Duv > nohup_06OnOff_Finding.out 2>&1 &
```
Example:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 6 UM_TKO_CIT_Duv --Indep > nohup_06OnOff_Finding_Indep350.out 2>&1 &
```

***Outputs:***  
Files are saved in `OnOff_target_ranking/` subdirectories created inside the input folder:
- ON-target file are saved as:  
    - `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_ON_target_<T><MC>.tsv`
    - `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_ON_target_<T><MC>_blackfilt.tsv`
- OFF-target file are saved as:
    - `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_OFF_target_<T><MC>.tsv`
    - `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspots_<MODE>_OFF_target_<T><MC>_blackfilt.tsv`

___

### **7) Top 200 Hotspots Selection** 
Option: `7` or `Top200_Selection`  

This step runs `07_top200hotspots_selection.sh`. It extracts the most significant hotspots (Top 200) based on the ranking performed in Step 4. 

It applies two different extraction methods:  
- Positional Extraction: Selects the first 200 data rows from the Full and Filtered files. 
- Rank-based Filtering: Uses awk to extract all regions where the value in the "Rank" column is $\le 200$.

***Usage:***
```bash
# Syntax:
./00_DepIndep_launch.sh 7 <ANALYSIS_SAMPLE_NAME> [FLAGS: --Dep / --Indep] [--mut / --dens]
```
Batch mode (recommended)
```bash
# Full analysis for a specific sample:
nohup ./00_DepIndep_launch.sh 7 UM_TKO_CIT_Duv > nohup_07top200_selection.out 2>&1 &
```
Example:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 7 UM_TKO_CIT_Duv --Dep > nohup_07top200_selection_Dep350.out 2>&1 &
```

***Outputs:***  
TSV files containing the Top 200 regions, identifying them with suffixes like `Ranked200hotspots` are saved in `Full_hotspts/` or in `Filtered_hotspots/`:
- `Full_hotspts/<E>_<SAMPLE>_<TYPE>_<CL>_Ranked200hotspots_<MODE>_ON_target_<T><MC>.tsv`
- `Filtered_hotspots/<E>_<SAMPLE>_<TYPE>_<CL>_Ranked200hotspots_<MODE>_ON_target_<T><MC>_blackfilt.tsv`

Additional two files are created and saved in `Filtered_hotspots/Filterd_and_blacklist_top200` that show the amount of filtered and deleted fragile hotspts in present in the full_Ranked200hotspots file:
- `<E>_<SAMPLE>_<TYPE>_<CL>_Ranked_hotspts_<MODE>_<T><MC>_200blackfilt`
- `<E>_<SAMPLE>_<TYPE>_<CL>_Deleted_hotspts_<MODE>_<T><MC>_200blacklist`

___

### **8) Coverage Table Generation (Detect-seq vs HTGTS)**
Option: `8` or `Hotspots_coverage`

This step runs `08_HTGTS-Detectseq_coverage.R`. It calculates the statistical overlap and coverage between Detect-seq results and HTGTS references. It generates **summary tables** comparing 
- "Detect-seq vs HTGTS" (How many Detect hotspots hit targets?) 
- "HTGTS vs Detect-seq" (Sensitivity analysis).

***Usage:***
```bash
# Syntax:
./00_DepIndep_launch.sh 8 <ANALYSIS_SAMPLE_NAME> [FLAGS: --Dep / --Indep] [--mut / --dens]
```
Batch mode (recommended)
```bash
# Full analysis for a specific sample:
nohup ./00_DepIndep_launch.sh 8 UM_TKO_CIT_Duv > nohup_08HTGTS-Detectseq_coverage.out 2>&1 &
```
Example:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 8 UM_TKO_CIT_Duv --Dep > nohup_08HTGTS-Detectseq_coverage_Dep350.out 2>&1 &
```
***Outputs:***  
Summary TSV files for both the Full dataset and Top 200 selection, located in:

```
|__ Full_hotspots
|    |__Pie_chart
|            |__Detect-seqvsHTGTS
|                   |__ ...TargetCoverageSummary_Detect.tsv
|                   |__ ...TargetCoverageSummary200_Detect.tsv
|            |__HTGTSvsDetect-seq
|                   |__ ...TargetCoverageSummary_HTGTS.tsv
|                   |__ ...TargetCoverageSummary200_HTGTS.tsv
|                   |__ ...TargetCoverageSummaryDetailed_HTGTS.tsv
|                   |__ ...TargetCoverageSummary200Detailed_HTGTS.tsv
|
|__ Filtered_hotspots
     |__Pie_chart
             |__Detect-seqvsHTGTS
                    |__ ...TargetCoverageSummary_Detect_blackfilt.tsv
                    |__ ...TargetCoverageSummary200_Detectblackfilt.tsv
             |__HTGTSvsDetect-seq
                    |__ ...TargetCoverageSummary_HTGTS_blackfilt.tsv
                    |__ ...TargetCoverageSummary200_HTGTS_blackfilt.tsv
                    |__ ...TargetCoverageSummaryDetailed_HTGTS_blackfilt.tsv
                    |__ ...TargetCoverageSummary200Detailed_HTGTS_blackfilt.tsv
```
___

### **9) Pie Chart Generation**

Option: `9` or `Coverage_Piechart`

This step runs `09_Coverage_Piechart.R`. It takes the summary tables generated in `Step 8` and produces publication-quality pie charts using ggplot2. The charts visualize the distribution of ON-target, OFF-target, and unknown target coverage.

***Usage:***
```bash
# Syntax:
./00_DepIndep_launch.sh 9 <ANALYSIS_SAMPLE_NAME> [FLAGS: --Dep / --Indep] [--mut / --dens]
```
Batch mode (recommended)
```bash
# Full analysis for a specific sample:
nohup ./00_DepIndep_launch.sh 9 UM_TKO_CIT_Duv > nohup_09Coverage_Piechart.out 2>&1 &
```
Example:
```bash
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 9 UM_TKO_CIT_Duv --Dep > nohup_09Coverage_Piechart_Dep350.out 2>&1 &
```

***Outputs:***
PNG images 
- `DetectseqvsHTGTS/`
    - `...CoveragePie_Detect...png`
    - `...CoveragePie200_Detect...png`
- `HTGTSvsDetectseq/`
    - `...CoveragePie_HTGTS...png`
    - `...CoveragePie200_HTGTS...png`
    - `...DetailedCoveragePie_HTGTS...png`
    - `...DetailedCoveragePie200_HTGTS...png`
___