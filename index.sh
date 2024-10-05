newvg="/home/projects/MAAG/Magpie/Magpie/vgan_corrected/dep/vg/bin/vg"

#### HAPLOCART

$newvg minimizer -o vgan_corrected/share/vgan/hcfiles/graph.min -t 50 -d vgan_corrected/share/vgan/hcfiles/graph.dist \
                -g vgan_corrected/share/vgan/hcfiles/graph.gbwt -p vgan_corrected/share/vgan/hcfiles/graph.og -k 10 -w 2

$newvg rymer -o vgan_corrected/share/vgan/hcfiles/graph.ry -t 50 -d vgan_corrected/share/vgan/hcfiles/graph.dist \
               -g vgan_corrected/share/vgan/hcfiles/graph.gbwt -p vgan_corrected/share/vgan/hcfiles/graph.og -k 10 -w 2
