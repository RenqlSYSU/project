#!/bin/bash

cd /home/ys17-23/Extension2/renql/TRACK-1.5.2/
filname=$1
output=$2
season=$3 # 0 = monthly; 1 = seasonal
echo $filname
echo $output

if [ $season == 0 ];then
    nday=(31 28 31 30 31 30 31 31 30 31 30 31)
    #nmonth=(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
    nmonth=($(seq 1 1 12))
    frae=744
fi

if [ $season == 1 ];then
    nday=(90 92 92 91)
    nmonth=(DJF MAM JJA SON)
    frae=0
fi

if [ $season == -1 ];then
    nday=(365)
    nmonth=(all)
    frae=744
fi


echo ${nmonth[*]}

input=indat/STATS.latlng_1hr.in
sed -i "21s/.*/${filname}/" $input 
for nm in $(seq 1 1 ${#nday[*]});do #{1..${nm}};do
    if [ ! -f ${output}_${nmonth[$((nm-1))]}.nc ];then
        fras=$((frae+1))
        frae=$((fras+24*nday[$((nm-1))]-1))
        echo "month=${nmonth[$((nm-1))]}, frame_s=${fras}, frame_e=${frae}"

        sed -i "34s/.*/${fras}/" $input 
        sed -i "35s/.*/${frae}/" $input
        
        bin/track.linux < $input > ${output}_record
        mv outdat/stat_trs_scl.linux_1.nc ${output}_${nmonth[$((nm-1))]}.nc
    fi
done

