#for file in data/*; do
#    nice ./vgan haplocart -fq1 $file -t 3 -p -o outputs/$file --hcfiles /home/projects/mito_haplotype/vgan/share/hcfiles/
#done

# Run on full BAMS
for file in full/*.bam; do
    # Extract the base filename (without extension) of the bam file
    base_filename=$(basename ${bam_file} .bam)

    nice samtools bam2fq $file | /home/projects/mito_haplotype/vgan/bin/vgan haplocart -fq1 /dev/stdin -t -1;
                        done
