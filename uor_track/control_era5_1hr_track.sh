#!/bin/bash
#SBATCH --partition=cluster
#SBATCH --output=/home/users/qd201969/ERA5-1HR-lev/TMP/slurm-%j.out 
#SBATCH --error=/home/users/qd201969/ERA5-1HR-lev/TMP/slurm-%j.err
#SBATCH --job-name=ew
#SBATCH --time=10-00:00:00

lev=(850 500 250)
OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
prefix=$1 #ff # tr is original cyclone ; ff is the filtered cyclone
# fft is the file that has been converted to real time

pro=$2 #2  # -1=convert to real time ;0 = run track ; 1 = combine ; 
# 2 = statistics ; 3 = box filt 
# 4 = three level match ;5 = two level match ; 6 = combine match file
filt=$3 #1 #filt=1, then use filter track to match
suffix=$4 #use to draw statistics, tpfilter, added on 20220601
nl1=0
nl2=1
nl3=0

# filter area
#option=$4 #5 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)/pass no genesis(5)
#lats=$5 #25 #25
#latn=$6 #45 #45
#lonl=$7 #60 #60
#lonr=$8 #110 #80
#if [ $filt == 1 ]; then 
#    suffix=_${option}_${lats}${latn}-${lonl}${lonr}
#else
#    suffix=""
#fi

if [ $pro == -1 ];then
    cd /home/users/qd201969/TRACK-1.5.2
    pwd
    echo "convert to real time"
    timestep=1 #hour

    for ny in {1980..2020};do
        starttime=$((ny-1))120100
        for nl in {0..2};do
            filname1=${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${prefix}_trs_pos.${suffix}
            filname2=${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${prefix}t_trs_pos.${suffix}
            if [ ! -f "${filname2}" ]; then
                utils/bin/count ${filname1} 0 0 5 4 0 $starttime $timestep
                mv ${filname1}.new ${filname2}
            fi
        done
    done
fi

if [ $pro == 0 ];then
    nl=2
    cd ~/TRACK-1.5.2/
    pwd
    echo "run track run_era5_1hr_qiaoling.csh ${lev[$nl]}"
    for ny in {2007..2020};do
        ./run_era5_1hr_qiaoling.csh ${ny} ${lev[$nl]} > ~/ERA5-1HR-lev/record_${lev[$nl]}_${ny}
    done
fi

if [ $pro == 1 ];then
    cd $OUTDIR
    pwd
    echo "combine ${prefix}_trs_pos"
    for nl in {0..2};do
        echo 41 > combine.in_${lev[$nl]}
        echo 1 >> combine.in_${lev[$nl]}
        
        for ny in {1980..2020};do
            filname=${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${prefix}_trs_pos.${suffix}
            echo $filname >> combine.in_${lev[$nl]}
        
            if [ ! -f "$filname" ]; then
                echo "$filname does not exist"
                #rm -irf ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET
                #~/TRACK-1.5.2/run_era5_1hr_qiaoling.csh ${ny} ${lev[$nl]}
            else
                #awk '{if(NF==4 && ($2 > 400 || $3 > 100 || $4 > 1000)) print FILENAME, $0}' \
                #    ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${prefix_trs_pos}
                num=$(awk '{if(NF==4 && ($2 > 400 || $3 > 100 || $4 > 1000)) print FILENAME, $0}' $filname | awk 'END{print NR}')
                if [ ${num} -ge 2 ]; then 
                    echo "Need to rerun ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET, whose error line is ${num}"
                    #rm -irf ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET
                    #~/TRACK-1.5.2/run_era5_1hr_qiaoling.csh ${ny} ${lev[$nl]}
                fi
            fi
        done
        
        ~/TRACK-1.5.2/utils/bin/combine < combine.in_${lev[$nl]} > record_combine
        mv ./combined_tr_trs ./${prefix}add_${lev[$nl]}_1980-2020
        awk 'NR==4 {print FILENAME, $2}' ./${prefix}add_${lev[$nl]}_1980-2020 | tee -a ${OUTDIR}number
    done
fi

if [ $pro == 2 ];then
    echo "=========== statistics ================"
    cd ${OUTDIR} #match${suffix}
    path=$(pwd)
    if [ ! -d statistic ];then
        mkdir statistic
    fi

    np=0
    level=1 # 2=total cyclone level,1=middle range filt cyclone, 0=small range filt and match cyclone level
    for file in ${prefix}_*_1980-2020${suffix} ; do
    #for file in ${prefix}_*_match ; do
        filname="\/home\/users\/qd201969\/ERA5-1HR-lev\/${file}"
        output=${OUTDIR}statistic/${file}_ptstat
        #filname="\/home\/users\/qd201969\/ERA5-1HR-lev\/match${suffix}\/${file}"
        #output=${OUTDIR}match${suffix}/statistic/${file}_stat
        
        echo $file
        season=0 # 0 = monthly; 1 = seasonal; -1 all year
        if [ ! -f ${output}.nc ];then
            sh ~/uor_track/stat_1hr_dec_jan.sh ${filname} ${output} ${season}
            cdo -r -copy ${output}_[1-9].nc ${output}_1[0-2].nc ${output}.nc 
            rm ${output}_[1-9].nc 
            rm ${output}_1[0-2].nc
        fi
        #if [ $filt == 1 ]; then 
        #    python ~/uor_track/2109-draw_stat_con_monthly_xr_mp.py \
        #        ${file} ${output}.nc ${level} ${lev[$np]}${suffix} \
        #        1 ${lats} ${latn} ${lonl} ${lonr}
        #else
        #    python ~/uor_track/2109-draw_stat_con_monthly_xr_mp.py \
        #        ${file} ${output}.nc 2 ${lev[$np]}${suffix} 0
        #fi
        np=$((np+1))
    done
    #python ~/uor_track/2109-draw_stat_con_xr_mp.py \
    #    ${prefix} ${suffix} ${level} 0 ${lats} ${latn} ${lonl} ${lonr}
    python ~/uor_track/2109-draw_stat_con_seasonal_xr_mp.py \
        ${prefix} ${suffix} ${level} 0 ${lats} ${latn} ${lonl} ${lonr}
fi

if [ $pro == 3 ];then
    echo "=========== trs filt (${lats}-${latn}N, ${lonl}-${lonr}E) ================"
    #box [filname] [Lat1.] [Lat2.] [Lng1.] [Lng2.] 
    # [Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)] 
    # [Negate (1)] [Track Thresh] [Start Time, YYYYMMDDHH]

    cd ~/TRACK-1.5.2/
    pwd
    for nl in {0..2};do
        file=${prefix}_${lev[$nl]}_1980-2020
        echo "lev: ${lev[$nl]} box filter"

        #python ~/uor_track/2110-line_filter.py ${OUTDIR}$file ${lats} ${latn} ${lonl} ${lonr} $option 24 
        if [ $option -le 4 ];then 
            utils/bin/box ${OUTDIR}$file ${lats} ${latn} ${lonl} ${lonr} $option 0 0.0 << EOF
n
EOF
            mv ${OUTDIR}${file}.new ${OUTDIR}${file}${suffix}
        else
            utils/bin/box ${OUTDIR}$file ${lats} ${latn} ${lonl} ${lonr} 2 0 0.0 << EOF
n
EOF
            mv ${OUTDIR}${file}.new ${OUTDIR}${file}${suffix}
            utils/bin/box ${OUTDIR}${file}${suffix} ${lats} ${latn} ${lonl} ${lonr} 0 1 0.0 <<EOF 
n
EOF
            mv ${OUTDIR}${file}${suffix}.new ${OUTDIR}${file}${suffix}
        fi
        awk 'NR==4 {print FILENAME, $2}' ${OUTDIR}${file}${suffix} | tee -a ${OUTDIR}number
        
#        utils/bin/tr2nc ${OUTDIR}${file}${suffix} s utils/TR2NC/tr2nc.meta
    done
#    cd ~/uor_track/fig/
#    python ~/uor_track/2107-draw_panel_traj.py $option ${lats} ${latn} ${lonl} ${lonr}
fi

if [ $pro == 4 ];then
    echo "=========== match low level with high level ==================="
    echo "Firts match 250 and 500, then match 500match250_yes and 850, match 500match250_no and 850"
    dist=5.0 # distance threshold
    timt=0.1 # time threshold

    cd ~/TRACK-1.5.2/
    pwd

    for ny in {1980..2020};do
    #dir=${OUTDIR}match${suffix}
    #file1=${OUTDIR}${prefix}_${lev[$nl1]}_1980-2020${suffix}
    #file2=${OUTDIR}${prefix}_${lev[$nl2]}_1980-2020${suffix}
    #file3=${OUTDIR}${prefix}_${lev[$nl3]}_1980-2020${suffix}
    
    dir=${OUTDIR}match_${ny}
    file1=${OUTDIR}ERA5_VOR${lev[$nl1]}_1hr_${ny}_DET/${prefix}_trs_pos
    file2=${OUTDIR}ERA5_VOR${lev[$nl2]}_1hr_${ny}_DET/${prefix}_trs_pos
    file3=${OUTDIR}ERA5_VOR${lev[$nl3]}_1hr_${ny}_DET/${prefix}_trs_pos

    if [ -d $dir ];then
        rm $dir/ff_*
    else
        mkdir $dir 
    fi

    utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s -1 1 ${dist} ${timt} > ${dir}/record${lev[$nl1]}_${lev[$nl2]}
    mv ./match_ens* $dir
    mv ./str.dat ${dir}/str.dat.${lev[$nl1]}_${lev[$nl2]}
    mv ./diff.dat ${dir}/diff.dat.${lev[$nl1]}_${lev[$nl2]}
    mv ./temp-dist.stat ${dir}/diff.dat.${lev[$nl1]}_${lev[$nl2]}
    rename match_ens1 ${prefix}_${lev[$nl1]}_${lev[$nl2]} $dir/*
    rename match_ens2 ${prefix}_${lev[$nl2]}_${lev[$nl1]} $dir/*
    
    utils/bin/censemble2 $file3 $file2 0 100 10 1 0 0 0 0 s -1 1 ${dist} ${timt} > ${dir}/record${lev[$nl3]}_${lev[$nl2]}
    mv ./match_ens* $dir
    mv ./str.dat ${dir}/str.dat.${lev[$nl3]}_${lev[$nl2]}
    mv ./diff.dat ${dir}/diff.dat.${lev[$nl3]}_${lev[$nl2]}
    mv ./temp-dist.stat ${dir}/diff.dat.${lev[$nl3]}_${lev[$nl2]}
    rename match_ens1 ${prefix}_${lev[$nl3]}_${lev[$nl2]} $dir/*
    rm ${dir}/match_ens2*
    
    term=(${lev[$nl2]}_${lev[$nl1]}_yes ${lev[$nl2]}_${lev[$nl1]}_no ${lev[$nl1]}_${lev[$nl2]}_yes)
    for pre in ${term[*]};do
        utils/bin/censemble2 ${dir}/${prefix}_${pre}.dat $file3 0 100 10 1 0 0 0 0 s -1 1 ${dist} ${timt} > ${dir}/record${pre}_${lev[$nl3]}
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
        awk 'NR==4 {print FILENAME, $2}' ${file} | tee -a ${OUTDIR}number
    done
    done
fi

if [ $pro == 5 ];then
    echo "=========== match low level with high level ==================="
    dist=5.0 # distance threshold
    timt=0.1 # time threshold
    cd ~/TRACK-1.5.2/
    pwd
    if [ $filt == 1 ]; then 
        #dir=${OUTDIR}${suffix}match${lev[$nl1]}_${lev[$nl2]}
        #file1=${OUTDIR}${prefix}_${lev[$nl1]}_1980-2020${suffix}
        #file2=${OUTDIR}${prefix}_${lev[$nl2]}_1980-2020${suffix}
        dir=${OUTDIR}pas_match${lev[$nl1]}_${lev[$nl2]}
        file1=/home/users/qd201969/ERA5-1HR-lev/ff_${lev[$nl1]}_1980-2020_2_2745-6060_PAS
        file2=/home/users/qd201969/ERA5-1HR-lev/ff_${lev[$nl2]}_1980-2020_2_2745-6060_PAS
    else
        dir=${OUTDIR}match${lev[$nl1]}_${lev[$nl2]}
        file1=${OUTDIR}${prefix}_${lev[$nl1]}_1980-2020
        file2=${OUTDIR}${prefix}_${lev[$nl2]}_1980-2020
    fi
    
    if [ -d $dir ];then
        rm $dir/*
    else
        mkdir $dir 
    fi

    #utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s 0 1 ${dist} ${timt} > ${dir}/record
    utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s -1 1 ${dist} ${timt} > ${dir}/record
    mv ./match_ens* $dir
    mv ./str.dat $dir
    mv ./diff.dat $dir
    mv ./temp-dist.stat $dir

    utils/bin/tr2nc ${dir}/match_ens1_yes.dat s utils/TR2NC/tr2nc.meta
    utils/bin/tr2nc ${dir}/match_ens1_no.dat s utils/TR2NC/tr2nc.meta
    utils/bin/tr2nc ${dir}/match_ens2_yes.dat s utils/TR2NC/tr2nc.meta
    utils/bin/tr2nc ${dir}/match_ens2_no.dat s utils/TR2NC/tr2nc.meta

    rename .dat.nc .nc $dir/*
    rename match_ens1 ${lev[$nl1]}_${lev[$nl2]} $dir/*
    rename match_ens2 ${lev[$nl2]}_${lev[$nl1]} $dir/*
fi

if [ $pro == 6 ];then
    echo "combine ${prefix}_trs_pos for yearly match results"
    cd ${OUTDIR}match_1980
    pwd
    for file in ${prefix}_*;do
        echo $file
        echo 41 > ${OUTDIR}combine.in_$file
        echo 1 >> ${OUTDIR}combine.in_$file
        
        for ny in {1980..2020};do
            echo ${OUTDIR}match_${ny}/${file} >> ${OUTDIR}combine.in_$file
        
            if [ ! -f "${OUTDIR}match_${ny}/${file}" ]; then
                echo "${OUTDIR}match_${ny}/${file} does not exist"
            else
                num=$(awk '{if(NF==4 && ($2 > 400 || $3 > 100 || $4 > 1000)) print FILENAME, $0}' ${OUTDIR}match_${ny}/${file} | awk 'END{print NR}')
                if [ ${num} -ge 2 ]; then 
                    echo "Need to rerun ${OUTDIR}match_${ny}/${file}, whose error line is ${num}"
                fi
            fi
        done
        
        ~/TRACK-1.5.2/utils/bin/combine < ${OUTDIR}combine.in_$file > ${OUTDIR}record
        mv ./combined_tr_trs ${OUTDIR}${file}_match
        awk 'NR==4 {print FILENAME, $2}' ${OUTDIR}${file}_match | tee -a ${OUTDIR}number
    done
fi
