# SAFARI Experiment Data

This repository contains the scripts, data, and analyses necessary to reproduce our benchmarking experiment results. Please note that you may need to modify some paths and configurations to reproduce the results on your machine.

## Scripts

- **Euka Experiment**: 
  - Script: `./run_euka.sh`
  - Description: Runs the euka experiment, analyzing eukaryotic sequence data.

- **HaploCart Experiment**: 
  - Script: `./run_hc.sh`
  - Description: Executes the HaploCart experiment, focusing on haplotype cartography.

- **Index Creation**: 
  - Script: `./index.sh`
  - Description: Generates the minimizer and RYmer index files on disk, critical for subsequent analyses.

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
  - Description: Produces simulations to benchmark linear aligners, underpinning comparative analyses.

- **Compute statistics from BAM file**:
  - Script: `parseBamMito.py`
  - Description: Extracts statistical data from BAM files, supporting comprehensive performance evaluation.

- **Damage Rate Estimation**:
  - Script: `damage.py`
  - Description: Generates a CSV file with results from the damage rate estimation experiment, offering insights into DNA degradation patterns.

- **Generate Substitution Matrices**:
  - Script: `make_subs.py`
  - Description: Creates substitution matrices per sample, facilitating the analysis of sequence variation and evolutionary dynamics.

## spurious_model Directory

- **Estimate parameters for spurious alignment model**:
  - Script: `main.py`
  - Description: Estimates parameters for the power law distribution of our spurious alignment model using RYmer and minimizer indexes.

## Additional Folders

- **Index Logs** (`index_logs`):
  - Contains detailed logs and data regarding the construction of minimizer and RYmer indexes with various (k,w) parameter sets, facilitating replication and modification of our indexing strategy.

- **Statistical Analysis and Plotting** (`linear_experiment/human_mito/stat_results`):
  - Houses scripts for statistical analysis and visualization of data from linear alignment and damage rate estimation experiments, enabling in-depth exploration and presentation of findings.

- **Alignment and Performance Data** (`linear_experiment/human_mito/alignments`):
  - Stores alignment outputs, performance statistics, and substitution matrices, serving as a repository for raw and processed data crucial to our study's insights.

