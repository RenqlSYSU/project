#!/bin/bash
#value=`cat ./data/JET_index_MEJindx_2040N5080E_7918_DJF.txt`
#echo $value[1]
#str1="(/"${value[*]}"/)"
#echo ${str1//${IFS:0:1}/,}

str1="(/"
first=`cat ./data/JET_index_MEJindx_2040N5080E_7918_DJF.txt | head -n 1`
str1=$str1$first
for line in `cat ./data/JET_index_MEJindx_2040N5080E_7918_DJF.txt | tail -n +2`
do
    str1=$str1","$line
done
str1=$str1"/)"
echo $str1
