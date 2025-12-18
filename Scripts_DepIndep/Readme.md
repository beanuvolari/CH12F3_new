# Define AID dependent and independent hotspots

## 00_DepIndep_launch

This repository contains a pipeline to define and visualize AID-dependent and independent hotspots. The main script ```00_DepIndep_launch.sh``` acts as a launcher to orchestrate the different steps of the analysis using Docker containers. This ensures modularity and reproducibility across different environments.

***Prerequisites***

- Docker must be installed and running.
- The Docker image detectseqpipe:4 (or the version specified in the script).

***Configuration***

To run this script, update the following variables at the top of the script: 

- BASE_DIR: the **absolute paths** to the project folder;
```
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
```
# Cell lines to process
CELL_LINE="cell_line1"
THRESHOLD="threshold"
MIN_COUNT="min_count"
ENLARGEMENT="enlargement"
```
- ORDERED_SAMPLES: array with the sample names to process, in the order in which they should be analyzed.
```
# Samples'order
ORDERED_SAMPLES=("sample_1" "sample_2" "sample_3" ...)
```


***Usage***  
```
# Syntax
./00_DepIndep_launch.sh <step> [options...]

# Help
./00_DepIndep_launch.sh --help
./00_DepIndep_launch.sh -h
```
- If no 'step' is specified, the script will prompt you an error. Have a look of the help.  
- The script automatically handles volume mapping between your host directories and the Docker container (/scratch/...)


## Steps


#### **1) Generation of a cross-sample shared hotspots dataset**

Option: ```1``` or ```Shared_Hotspots```

This step runs ```01_0SharedHotspots.sh```. It wraps ```01_1SharedHotspots.R``` that collects and merges hotspot regions ("RIDER" regions) detected across multiple samples, producing a consolidated dataset. It allows for the optional expansion of hotspot coordinates via the **enlargement** parameter.


***Input Requirement:*** The script expects the RIDER output:
- Location: ```/scratch/pmat/Filtered_pmat/<SAMPLE>/merged_CTGA/threshold_<T>/<T>_rider<MC>```  
-- Filename format: ```<sample_name>_merged_CTGA_<T>_sorted-RIDER.clean*.bed```


***Usage:***
```
# Run:
./00_DepIndep_launch.sh 1

# Run via command line (parameters are defined above in the script config part)

# Syntax:
nohup ./00_DepIndep_launch.sh 1 > nohup_01SharedHotspots_merge.out 2>&1 &
```
Example:
```
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh 1 > nohup_01SharedHotspots_merge350.out 2>&1 &
```


***Output:***  
The output table is saved in 
- Location: ```DepIndep_dataset/Threshold_<T>/Enlargement_<E>/```

- Filename format: ```<E>_<CELL>_SharedHotspots_<T>_rider<MC>.bed``` 

The script is designed to be run in a **batch mode**.  

___

#### **2) Plotting (Euler, Upset, Venn)**

Option: ```2``` or ```Plots```

This step runs ```02_0SharedHotspots_plots.sh``` (which wraps three R scripts: ```02_Euler_plot.R```, ```02_Upset_plot.R``` and ```02_Venn_plot.R```). It generates visualizations to analyze the intersections of the hotspots defined in Step 1.  

***Input Requirement:*** The script expects a tab-delimited BED-like file generated in Step 1.

- Location: ```DepIndep_dataset/Threshold_<T>/Enlargement_<E>/```

- Filename format: ```<E>_<CELL>_SharedHotspots_<T>_rider<MC>.bed```  


***Usage & Flags:***

You can run all plots at once or select specific ones using flags.

```
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
```
# Run via command line (parameters are defined above in the script config prt)

# Generate ALL plots (Default):
nohup ./00_DepIndep_launch.sh > nohup_DepIndep_plot.out 2>&1 &

# Generate only Euler and UpSet plots:
nohup ./00_DepIndep_launch.sh 2 -e -u > nohup_DepIndep_ploteu.out 2>&1 &
```

Examples:
```
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0

nohup ./00_DepIndep_launch.sh > nohup_DepIndep_plot350.out 2>&1 &

nohup ./00_DepIndep_launch.sh 2 -e -u > nohup_DepIndep_ploteu350.out 2>&1 &
```


If you run ```./00_DepIndep_launch.sh 2``` without flags, the script will generate all the plot types.

***Outputs:*** Plots are saved as high-resolution PNGs (300 DPI) in: ```DepIndep_dataset/Threshold_<T>/Enlargement_<E>/Plots_Euler_Venn_Upset/```
- Euler: ```<E>_Plots_Euler_<T>_<MC>.png```
- Venn:  ```<E>_Plots_Venn_<T>_<MC>.png```
- Upset: ```<E>_Plots_Upset_<T>_<MC>.png```


***Note on R Scripts:***

- Venn: will automatically skip generation if the input dataset contains more than 5 sample columns.

- Euler: uses the eulerr package. It may fail if the set intersections are too complex to map geometrically.

- UpSet: converts "+" values to binary (0/1) for intersection analysis.  

***Troubleshooting***  
- ```"[ERROR] Your input file shows formatting problems"```: The plotting script checks if the input file is strictly tab-delimited. If this error occurs, ensure your dataset does not contain spaces instead of tabs or malformed line endings (e.g., Windows \r).

- Venn Plot missing: Check the logs. If your data has >5 sample columns, the Venn diagram is intentionally skipped.
___

#### **3) AID-Dependent/ AID-Independent hotspots split**

Option: ```3 ```or ```AID_DepIndep_Split``` 

This step runs the script ```03_AID_DepIndephotspots_split.sh```, which separates the shared hotspot dataset (generated in Step 1) into two categories for a user-selected sample:  
- **AID-Dependent** hotspots: regions where the selected sample carries a “**–**”

- **AID-Independent** hotspots: regions where the selected sample carries a “**+**”

The script identifies the correct sample column in the input file and automatically generates two BED-like outputs.

**Input Requirement:**
This step expects as input the shared hotspots file created in Step 1:

Location:
```DepIndep_dataset/Threshold_<T>/Enlargement_<E>/```

Filename format:
```<E>_<CELL>_SharedHotspots_<T>_rider<MC>.bed```


***Parameters:***

You must supply one parameter: 
- REF_SAMPLE_NAME – The specific sample whose “+/–” column will define Dependent vs Independent hotspots.


***Usage:***
```
# Syntax:
./00_DepIndep_launch.sh 3 <REF_SAMPLE_NAME>
```

Batch mode (recommended)
```
# Run via command line (parameters are defined above in the script config part)

# Syntax:
nohup ./00_DepIndep_launch.sh 3 <REF_SAMPLE_NAME> > nohup_03AID_DepIndephotspots_split.out 2>&1 &
```

Examples:
```
# cell_line: CH12F3, theshold:3, min_count:5, enlargement:0
nohup ./00_DepIndep_launch.sh AID_KO > nohup_03AID_DepIndephotspots_split350.out 2>&1 &
```

***Outputs:*** Two table are generate
- Dependent:```/scratch/DepIndep_dataset/Threshold_<T>/Enlargement_<E>/Dependent_hotspots/<E>_<CL>_AID_Dependent_hotspots_<T>_<MC>.bed```
- Independent:```/scratch/DepIndep_dataset/Threshold_<T>/Enlargement_<E>/Independent_hotspots/<E>_<CL>_AID_Independent_hotspots_<T>_<MC>.bed```
___
