#!/bin/bash

#lev=(850 500 250)
OUTDIR=/home/ys17-23/Extension2/renql/project/uor_track/mdata/
files[1]=ff_match_850local_500remote_6dist
files[2]=ff_match_500local_250remote_6dist
files[3]=ff_match_250local_500remote_6dist
files[4]=ff_match_850local_500local_6dist
files[5]=ff_match_500local_250local_6dist
files[6]=ff_match_250local_500local_6dist
outfiles[1]=ff_match_850localremote_6dist
outfiles[2]=ff_match_500localremote_6dist
outfiles[3]=ff_match_250localremote_6dist
outfiles[4]=ff_match_850locallocal_6dist
outfiles[5]=ff_match_500locallocal_6dist
outfiles[6]=ff_match_250locallocal_6dist

echo "=========== statistics ================"
cd ${OUTDIR} #match${suffix}
path=$(pwd)
if [ ! -d statistic ];then
    mkdir statistic
fi

np=0
level=0 # 2=total cyclone level,1=middle range filt cyclone, 0=small range filt and match cyclone level
for nf in {1..6} ; do
    filname="${OUTDIR}${files[nf]}"
    output=${OUTDIR}statistic/${outfiles[nf]}
    
    echo $filname
    season=0 # 0 = monthly; 1 = seasonal; -1 all year
    if [ ! -f ${output}.nc ];then
        bash /home/ys17-23/Extension2/renql/project/uor_track/stat_1hr_dec_jan.sh ${filname} ${output} ${season}
        cdo -r -copy ${output}_[1-9].nc ${output}_1[0-2].nc ${output}.nc 
        rm ${output}_[1-9].nc 
        rm ${output}_1[0-2].nc
    fi
done
#python /home/ys17-23/Extension2/renql/project/uor_track/2109-draw_stat_con_xr_mp.py \
#    'ff' 'local' ${level} 0 0 0 0 0 
python /home/ys17-23/Extension2/renql/project/uor_track/2109-draw_stat_con_seasonal_xr_mp.py \
    'ff' 'local' ${level} 0 0 0 0 0 
#python /home/ys17-23/Extension2/renql/project/uor_track/2109-draw_stat_con_xr_mp.py \
#    'ff' 'remote' ${level} 0 0 0 0 0 
python /home/ys17-23/Extension2/renql/project/uor_track/2109-draw_stat_con_seasonal_xr_mp.py \
    'ff' 'remote' ${level} 0 0 0 0 0 

