# SAFARI Experiment Data

This repository contains the scripts, data, and analyses necessary to reproduce our benchmarking experiment results. Please note that you may need to modify some paths and configurations to reproduce the results on your machine.

## Scripts

- **Euka Experiment**: 
  - Script: `./run_euka.sh`
  - Description: Runs the euka experiment.

- **HaploCart Experiment**: 
  - Script: `./run_hc.sh`
  - Description: Runs the HaploCart experiment.

- **Index Creation**: 
  - Script: `./index.sh`
  - Description: Generates the minimizer and RYmer index files on disk.

## linear_experiment/human_mito Directory

- **Evaluate**: 
  - Script: `evaluate.py`
  - Description: Computes summary statistics from the `alignment_stats.csv` file for alignment quality assessment.

- **Compare Damage Estimate**: 
  - Script: `compare_damage_est.py`
  - Description: Performs cross-comparison of damage rate estimates derived from bam2prof analyses.

- **Interleave FASTQ**: 
  - Script: `interleave_fastq.sh`
  - Description: Interleaves FASTQ files, preparing them for alignment.

- **Generate simulations for linear aligner benchmarking experiment**:
  - Script: `Snakefile`
  - Description: Produces simulations to benchmark linear aligners compared to giraffe and SAFARI.

- **Compute statistics from BAM file**:
  - Script: `parseBamMito.py`
  - Description: Extracts alignment stats from BAM files.

- **Damage Rate Estimation**:
  - Script: `damage.py`
  - Description: Generates a CSV file with results from the damage rate estimation experiment.

## spurious_model Directory

- **Estimate parameters for spurious alignment model**:
  - Script: `main.py`
  - Description: Estimates parameters for the power law distribution of our spurious alignment model using RYmer and minimizer indexes.

## Additional Folders

- **Index Logs** (`index_logs`):
  - Contains detailed logs and data regarding the construction of minimizer and RYmer indexes with various (k,w) parameter sets.

- **Statistical Analysis and Plotting** (`linear_experiment/human_mito/stat_results`):
  - Houses scripts for statistical analysis and visualization of data from linear alignment and damage rate estimation experiments.

- **Alignment and Performance Data** (`linear_experiment/human_mito/alignments`):
  - Stores alignment outputs, performance statistics, and substitution matrices.

