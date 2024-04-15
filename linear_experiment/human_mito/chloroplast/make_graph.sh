vg_path="/home/projects/MAAG/Magpie/Magpie/vgan_corrected/dep/vg/bin/vg"

$vg_path construct -p -m 8 -M aribidopsis.fa > chloroplast.vg
$vg_path convert -t 60 -o chloroplast.vg > chloroplast.og
$vg_path view chloroplast.og > chloroplast.gfa
$vg_path snarls chloroplast.og > chloroplast.snarls
$vg_path index -j chloroplast.dist chloroplast.og
$vg_path index -x chloroplast.xg chloroplast.og
$vg_path gbwt -g chloroplast.giraffe.gbz --gbz-format -G chloroplast.gfa
$vg_path gbwt -o chloroplast.gbwt -g chloroplast.gg -Z chloroplast.giraffe.gbz
$vg_path minimizer -o chloroplast.min -g chloroplast.gbwt -d chloroplast.dist chloroplast.og
$vg_path rymer -o chloroplast.ry -g chloroplast.gbwt -d chloroplast.dist chloroplast.og
