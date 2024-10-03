
#./vgan_corrected/dep/vg/bin/vg minimizer -t 40 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist \
#-o ../euka_dir/euka_db.min ../euka_dir/euka_db.og -k 30 -w 20 && echo -e "\n\n"
#./vgan_corrected/dep/vg/bin/vg rymer -t 40 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist \
#-o ../euka_dir/euka_db.ry ../euka_dir/euka_db.og -k 30 -w 20 && echo -e "\n\n"

parallel -j 1 'zcat /home/projects2/euka_environments/full_cave_dhigh.fq.gz | /usr/bin/time -v ./vgan_corrected/bin/vgan euka -t 1 \
               -fq1 /dev/stdin --euka_dir /home/projects/MAAG/Magpie/euka_dir/ --minFrag 100 --minBins 6 \
               --deam3p /home/projects/MAAG/Magpie/Magpie/dhigh3p.prof --deam5p /home/projects/MAAG/Magpie/Magpie/dhigh5p.prof \
              -o euka_corrected_full_j_{} -j {}' ::: $(seq 0.9 0.1 0.9)

#zcat /home/projects2/euka_environments/full_cave_dhigh.fq.gz | ./vgan_uncorrected/bin/vgan euka -t 40 --minFrag 100 \
#-fq1 /dev/stdin --euka_dir /home/projects/MAAG/Magpie/euka_dir/ -o euka_uncorrected_full --minBins 6 --minMQ 10 --entropy 0.5


