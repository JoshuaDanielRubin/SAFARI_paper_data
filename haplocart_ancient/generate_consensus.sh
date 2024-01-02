
ANGSD_PATH=/home/ctools/angsd-0.935/angsd

for bam_file in full/*.bam; do

    # Extract the base filename (without extension) of the bam file
    base_filename=$(basename ${bam_file} .bam)

    samtools view -b full/$base_filename.bam MT -@ 10 > consensus/$base_filename.bam
    $ANGSD_PATH -seed 42 -doFasta 2 -minq 25 -minmapq 25 -uniqueonly 1 -out consensus/$base_filename -i consensus/$base_filename.bam -doCounts 1

    done
