for bam_file in full/*.bam; do

    # Extract the base filename (without extension) of the bam file
    base_filename=$(basename ${bam_file} .bam)

    ./haplocheck/cloudgene run haplocheck@1.3.2 --files /home/projects/MAAG/Josh/haplocart_ancient/full/individual_bams/$base_filename --format bam --threads 60 --output /home/projects/MAAG/Josh/haplocart_ancient/$base_filename

    done;
