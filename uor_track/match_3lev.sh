#!/bin/bash
prefix=$1 #ff or tr
suffix=$2 #_2_2745-6060_LYS

OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
lev=(250 500 850)
nl1=0
nl2=1
nl3=2
dist=5.0 # distance threshold
timt=0.1 # time threshold
reg_type=0 
#Region matching type, '-1' (no region matching), '0' (track in region),
#'1' (extrema in region), '2' (60% of points in region),'3' (extrema of region)

echo "=========== match low level with high level ==================="
cd ~/TRACK-1.5.2/
pwd

dir=${OUTDIR}match${suffix}
file1=${OUTDIR}${prefix}_${lev[$nl1]}_1980-2020${suffix}
file2=${OUTDIR}${prefix}_${lev[$nl2]}_1980-2020${suffix}
file3=${OUTDIR}${prefix}_${lev[$nl3]}_1980-2020${suffix}

if [ -d $dir ];then
    rm $dir/number
else
    mkdir $dir 
fi

utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s ${reg_type} 1 ${dist} ${timt} > ${dir}/record${lev[$nl1]}_${lev[$nl2]}
mv ./match_ens* $dir
mv ./str.dat ${dir}/str.dat.${lev[$nl1]}_${lev[$nl2]}
mv ./diff.dat ${dir}/diff.dat.${lev[$nl1]}_${lev[$nl2]}
mv ./temp-dist.stat ${dir}/diff.dat.${lev[$nl1]}_${lev[$nl2]}
rename match_ens1 ${prefix}_${lev[$nl1]}_${lev[$nl2]} $dir/*
rename match_ens2 ${prefix}_${lev[$nl2]}_${lev[$nl1]} $dir/*

utils/bin/censemble2 $file3 $file2 0 100 10 1 0 0 0 0 s ${reg_type} 1 ${dist} ${timt} > ${dir}/record${lev[$nl3]}_${lev[$nl2]}
mv ./match_ens* $dir
mv ./str.dat ${dir}/str.dat.${lev[$nl3]}_${lev[$nl2]}
mv ./diff.dat ${dir}/diff.dat.${lev[$nl3]}_${lev[$nl2]}
mv ./temp-dist.stat ${dir}/diff.dat.${lev[$nl3]}_${lev[$nl2]}
rename match_ens1 ${prefix}_${lev[$nl3]}_${lev[$nl2]} $dir/*
rm ${dir}/match_ens2*

term=(${lev[$nl2]}_${lev[$nl1]}_yes ${lev[$nl2]}_${lev[$nl1]}_no ${lev[$nl1]}_${lev[$nl2]}_yes)
for pre in ${term[*]};do
    utils/bin/censemble2 ${dir}/${prefix}_${pre}.dat $file3 0 100 10 1 0 0 0 0 s ${reg_type} 1 ${dist} ${timt} > ${dir}/record${pre}_${lev[$nl3]}
    mv ./match_ens* $dir
    mv ./str.dat ${dir}/str.dat.${pre}_${lev[$nl3]}
    mv ./diff.dat ${dir}/diff.dat.${pre}_${lev[$nl3]}
    mv ./temp-dist.stat ${dir}/diff.dat.${pre}_${lev[$nl3]}
    rename match_ens1 ${prefix}_${pre}_${lev[$nl3]} $dir/*
    rm ${dir}/match_ens2*
    rm ${dir}/${prefix}_${pre}.dat
done

pre=${lev[$nl3]}_${lev[$nl2]}_yes
utils/bin/censemble2 ${dir}/${prefix}_${pre}.dat $file1 0 100 10 1 0 0 0 0 s -1 1 ${dist} ${timt} > ${dir}/record${pre}_${lev[$nl1]}
mv ./match_ens* $dir
mv ./str.dat ${dir}/str.dat.${pre}_${lev[$nl1]}
mv ./diff.dat ${dir}/diff.dat.${pre}_${lev[$nl1]}
mv ./temp-dist.stat ${dir}/diff.dat.${pre}_${lev[$nl1]}
rename match_ens1 ${prefix}_${pre}_${lev[$nl1]} $dir/*
rm ${dir}/match_ens2*
rm ${dir}/${prefix}_${pre}.dat

rename ".dat" "" ${dir}/${prefix}_*

for file in ${dir}/${prefix}_*;do
    awk 'NR==4 {print FILENAME, $2}' ${file} | tee -a ${dir}/number
done

