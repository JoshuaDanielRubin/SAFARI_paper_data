/home/ctools/gargammel/art_src_MountRainier/art_illumina -i rCRS.fa -ss HS25 -na -o input -l 100 -f 1000

safari_path="/home/projects/MAAG/Magpie/Magpie/vgan_corrected/dep/vg/bin/vg"
giraffe_path="/home/projects/MAAG/Magpie/Magpie/vgan_uncorrected/dep/vg/bin/vg"

#/usr/bin/time -v $giraffe_path giraffe -p -m graph.min -t 60 -f input.fq -x graph.og -Z graph.giraffe.gbz -d graph.dist > giraffe.gam

/usr/bin/time -v $safari_path giraffe \
--deam-3p ../damageProfiles/none.prof \
--deam-5p ../damageProfiles/none.prof \
-p -m graph.min -q graph.ry -t 60 -f input.fq -x graph.og -Z graph.giraffe.gbz -d graph.dist > safari.gam
