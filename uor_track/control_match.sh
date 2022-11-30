#!/bin/bash

radiu=(3 9 12)
cd /home/ys17-23/Extension2/renql/project/uor_track/

for ra in ${radiu[@]} ; do
    python ./2210-match_local_local.py ${ra}
    python ./2210-match_local_remote.py ${ra}
    #python ./2211-convert_fftadd2ffadd.py ${ra}
    #bash ./control_stat_match.sh ${ra}
    #python ./2207-draw_stat_seasonal_shad_contour.py ${ra}
done
