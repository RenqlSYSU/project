#!/bin/bash
lev=(850 500 300)

cd ~/TRACK-1.5.2/
# filter area
lats=25
latn=50
lonl=40
lonr=70
option=0 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)

pro=2
filt=0 #if pro=2, filt=1, then use filter track to match
nl1=2
nl2=1

if [ $pro == 0 ];then
export PATH=${PATH}:.
pwd
echo $PATH
for nl in {2..2};do
    echo "lev: ${lev[$nl]}"

    sed -i "8s/.*/${lev[$nl]}00/" specfilt_T42.in

    bin/track.linux -i ERA5_VOR_6hr_2000_DET.nc -f file < specfilt_T42.in
    mv outdat/specfil.file_band001 indat/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42filt.dat
    rm outdat/specfil.file_band000

    master -c=ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42 -d=now -e=track.linux \
        -i=ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42filt.dat -f=f2000 -j=RUN_AT.in \
        -k=initial.T42_NH -n=1,64,24 -o=~/TRACK_result -r=RUN_AT_ \
        -s=RUNDATIN.VOR > record_master
    
    gunzip ~/TRACK_result/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42/tr_trs_pos.gz
    utils/bin/tr2nc ~/TRACK_result/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42/tr_trs_pos s utils/TR2NC/tr2nc.meta
    
    gunzip ~/TRACK_result/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42/ff_trs_pos.gz
    utils/bin/tr2nc ~/TRACK_result/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42/ff_trs_pos s utils/TR2NC/tr2nc.meta
done
fi

if [ $pro == 1 ];then
echo "=========== trs filt (${lats}-${latn}N, ${lonl}-${lonr}E) ================"
for nl in {0..2};do
    echo "lev: ${lev[$nl]}"
    file=~/TRACK_result/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42/tr_trs_pos

    utils/bin/box $file ${lats} ${latn} ${lonl} ${lonr} $option 0 0.0 
    mv ${file}.new ${file}.${option}_${lats}${latn}-${lonl}${lonr}
    utils/bin/tr2nc ${file}.${option}_${lats}${latn}-${lonl}${lonr} s utils/TR2NC/tr2nc.meta
done
fi

if [ $pro == 2 ];then
    echo "=========== match low level with high level ==================="
    if [ $filt == 1 ]; then 
        dir=~/TRACK_result/filt${option}_${lats}${latn}-${lonl}${lonr}match${lev[$nl1]}_${lev[$nl2]}
        file1=~/TRACK_result/ERA5_VOR${lev[$nl1]}_6hr_2000_DET_T42/tr_trs_pos.${option}_${lats}${latn}-${lonl}${lonr}
        file2=~/TRACK_result/ERA5_VOR${lev[$nl2]}_6hr_2000_DET_T42/tr_trs_pos.${option}_${lats}${latn}-${lonl}${lonr}
    else
        dir=~/TRACK_result/test_match${lev[$nl1]}_${lev[$nl2]}
        file1=~/TRACK_result/ERA5_VOR${lev[$nl1]}_6hr_2000_DET_T42/tr_trs_pos
        file2=~/TRACK_result/ERA5_VOR${lev[$nl2]}_6hr_2000_DET_T42/tr_trs_pos
    fi
    mkdir $dir 
    utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s -1 1 5.0 0.1 > ${dir}/record
    
    mv ./match_ens* $dir
    mv ./str.dat $dir
    mv ./diff.dat $dir
    mv ./temp-dist.stat $dir

    utils/bin/tr2nc ${dir}/match_ens1_yes.dat s utils/TR2NC/tr2nc.meta
    utils/bin/tr2nc ${dir}/match_ens1_no.dat s utils/TR2NC/tr2nc.meta
    utils/bin/tr2nc ${dir}/match_ens2_yes.dat s utils/TR2NC/tr2nc.meta
    utils/bin/tr2nc ${dir}/match_ens2_no.dat s utils/TR2NC/tr2nc.meta

    rename .dat.nc .nc $dir/*
    rename ens1 ${lev[$nl1]} $dir/*
    rename ens2 ${lev[$nl2]} $dir/*
fi

