
safari_path="/home/projects/MAAG/Magpie/Magpie/vgan_corrected/dep/vg/bin/vg"
giraffe_path="/home/projects/MAAG/Magpie/Magpie/vgan_uncorrected/dep/vg/bin/vg"
samtools_path="/home/ctools/bin/samtools"

for RATE in 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00; do
$safari_path giraffe -t 30 -f sample.fq.gz \
           -o BAM \
           -m ../../../rCRS_graph/graph.min \
           -d ../../../rCRS_graph/graph.dist \
           -q ../../../rCRS_graph/graph.ry \
           -Z ../../../rCRS_graph/graph.giraffe.gbz \
           --ref-path rCRS_PATHS \
           --deam-3p ATP2_mt_${RATE}.3p.prof --deam-5p ATP2_mt_${RATE}.5p.prof | $samtools_path sort > safari_rCRS_${RATE}.bam
done

for RATE in 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00; do
$safari_path giraffe -t 30 -f sample.fq.gz \
           -o BAM \
           -m /home/projects2/hominin/vgan_dev/share/tmfiles/graph.min \
           -d /home/projects2/hominin/vgan_dev/share/tmfiles/graph.dist \
           -q /home/projects2/hominin/vgan_dev/share/tmfiles/graph.ry \
           -Z /home/projects2/hominin/vgan_dev/share/tmfiles/graph.giraffe.gbz \
           --ref-path hominin_PATHS \
           --deam-3p ATP2_mt_${RATE}.3p.prof --deam-5p ATP2_mt_${RATE}.5p.prof | $samtools_path sort > safari_hominin_${RATE}.bam
done

