#!/bin/bash
echo "=========== match low level with high level ==================="
cd ~/TRACK-1.5.2/
dist=5.0 # distance threshold
timt=0.1 # time threshold
dir=$(pwd)
file1=/home/users/qd201969/ERA5-1HR-lev/ff_${lev[$nl1]}_1980-2020_2_2745-6060_PAS
file2=/home/users/qd201969/ERA5-1HR-lev/ff_${lev[$nl2]}_1980-2020_2_2745-6060_PAS

utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s 0 1 ${dist} ${timt} > ${dir}/record

utils/bin/tr2nc ${dir}/match_ens1_yes.dat s utils/TR2NC/tr2nc.meta
utils/bin/tr2nc ${dir}/match_ens1_no.dat s utils/TR2NC/tr2nc.meta
utils/bin/tr2nc ${dir}/match_ens2_yes.dat s utils/TR2NC/tr2nc.meta
utils/bin/tr2nc ${dir}/match_ens2_no.dat s utils/TR2NC/tr2nc.meta

rename .dat.nc .nc $dir/*
rename match_ens1 ${lev[$nl1]}_${lev[$nl2]} $dir/*
rename match_ens2 ${lev[$nl2]}_${lev[$nl1]} $dir/*

