#!/bin/bash

lev=(850 500 250)
OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
prefix=ff # tr is original cyclone ; ff is the filtered cyclone

filt=0 #filt=1, then use filter track to match

# filter area
lats=25
latn=45
lonl=60
lonr=80
option=2 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
if [ $filt == 1 ]; then 
    suffix=_${option}_${lats}${latn}-${lonl}${lonr}
else
    suffix=
fi

#figtitle=('only_250' '250-500' '250-500-850' \
#         'only_500' '500-850' '500-250' '500-250-850' \
#         'only_850' '850-500' '850-500-250')
figtitle=('250' '500' '850')

cd ${OUTDIR}
#cd ${OUTDIR}match${suffix}/
level=2 # 2=total cyclone level,1=middle range filt cyclone, 0=small range filt and match cyclone level
path=$(pwd)

np=0
for filename in ${prefix}_*_1980-2020${suffix};do
#for filename in ${prefix}_*_match;do
#for filename in ${prefix}_*;do # mainly used in ${OUTDIR}match${suffix}/
    echo ${filename}
    file=${path}/statistic/${filename}_stat

    #python ~/uor_track/2108-draw_stat_con_monthly.py \
    python ~/uor_track/2109-draw_stat_con_monthly_xr_mp.py \
        ${filename} ${file}.nc ${level} ${figtitle[$np]}${suffix} \
        0 ${lats} ${latn} ${lonl} ${lonr}
    np=$((np+1))
done

#mogrify -bordercolor white -trim ./month_*.png

