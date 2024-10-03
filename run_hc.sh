#!/bin/bash

# Specify the total number of threads
total_threads=1

# Specify the list of subsampling rates
#rates=("0.25x" "0.5x" "1x", "2x")
rates=("0.25x")

process_file_corrected() {
    bam_file=$1
    threads_per_job=$2
    rate=$3

    # Get the base name of the bam file (without extension)
    base_name=$(basename "$bam_file" .bam)

    # Define the expected output file path (modify as needed)
    output_file="hc_results/$base_name.corrected.log"

    # Check if the output file already exists
    if [ -f "$output_file" ]; then
        echo "Output file $output_file already exists, skipping $bam_file" >&2
        return
    fi

    # Log the initiation of the process for this bam file to stderr
    echo "Processing $bam_file" >&2

    # Run the haplocart pipeline
    samtools bam2fq "$bam_file" -@ 1 | /home/projects/MAAG/Magpie/Magpie/vgan_corrected/bin/vgan haplocart -np -t 1 -fq1 /dev/stdin \
    --hc-files /home/projects/MAAG/Magpie/Magpie/vgan_corrected/share/vgan/hcfiles \
    --deam3p /home/projects/MAAG/Magpie/Magpie/dhigh3p.prof \
    --deam5p /home/projects/MAAG/Magpie/Magpie/dhigh5p.prof \
    &>> "$output_file" 2>&1
}

process_file_uncorrected() {
    bam_file=$1
    threads_per_job=$2
    rate=$3

    # Get the base name of the bam file (without extension)
    base_name=$(basename "$bam_file" .bam)

    # Log the initiation of the process for this bam file to stderr
    echo "Processing $bam_file" >&2

    # Run the haplocart pipeline
    samtools bam2fq "$bam_file" -@ 20 | /home/projects/MAAG/Magpie/Magpie/vgan_uncorrected/bin/vgan haplocart -np -t 35 -fq1 /dev/stdin \
    --hc-files /home/projects/MAAG/Magpie/Magpie/vgan_corrected/share/vgan/hcfiles \
    &>> "hc_results/$base_name.uncorrected.log" 2>&1

}

# Export the function to be available in the parallel environment
export -f process_file_corrected
export -f process_file_uncorrected

# Loop over all subsampling rates
for rate in "${rates[@]}"; do
    # Specify the directory containing the bam files
    dir="/home/projects/MAAG/Magpie/Magpie/haplocart_ancient/subsampled_reps/$rate/"

    # Find the number of bam files
    num_files=$(ls $dir/*.bam | wc -l)

     # Calculate threads per job
    if [ $num_files -le $total_threads ]; then
        threads_per_job=1
    else
        threads_per_job=$(( total_threads / num_files ))
        threads_per_job=$(( threads_per_job < 1 ? 1 : threads_per_job ))
    fi

    # Run the function in parallel over all bam files
    nice -19 parallel --keep-order -j 1 process_file_corrected ::: $(ls $dir/*.bam) ::: $threads_per_job ::: $rate
    #nice -19 parallel --keep-order -j 16 process_file_uncorrected ::: $(ls $dir/*.bam) ::: $threads_per_job ::: $rate


done

