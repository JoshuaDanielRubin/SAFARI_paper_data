
#./vg_corrected/bin/vg minimizer -t 60 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist \
#-o ../euka_dir/euka_db.min ../euka_dir/euka_db.og -k 30 -w 10 && echo -e "\n\n"
#./vg_corrected/bin/vg rymer -t 60 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist \
#-o ../euka_dir/euka_db.ry ../euka_dir/euka_db.og -k 30 -w 10 && echo -e "\n\n"

#for POSTERIOR in 0.9;
#do
#    zcat /home/projects2/euka_environments/full_cave_dhigh.fq.gz | /usr/bin/time -v ./vgan_corrected/bin/vgan euka \
#    -t 20 -fq1 /dev/stdin --euka_dir /home/projects/MAAG/Magpie/euka_dir/ --minFrag 100 \
#    --deam3p /home/projects/MAAG/Magpie/Magpie/dhigh3p.prof --deam5p /home/projects/MAAG/Magpie/Magpie/dhigh5p.prof -o euka_corrected_full_j_${POSTERIOR} -j ${POSTERIOR}
#done

zcat /home/projects2/euka_environments/full_cave_dhigh.fq.gz | ./vgan_uncorrected/bin/vgan euka -t 20 --minFrag 100 \
-fq1 /dev/stdin --euka_dir /home/projects/MAAG/Magpie/euka_dir/ -o euka_uncorrected_full


### TESTING
#for POSTERIOR in 0.0;
#do
#./vgan_corrected/bin/vgan euka -t -1 -fq1 seq.fq --euka_dir /home/projects/MAAG/Magpie/euka_dir/ --deam3p /home/projects/MAAG/Magpie/Magpie/dhigh3p.prof --deam5p /home/projects/MAAG/Magpie/Magpie/dhigh5p.prof -o euka_corrected_full_j_${POSTERIOR} -j ${POSTERIOR}
#done
