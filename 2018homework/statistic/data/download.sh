#!/bin/bash
#====================================================================
# The purpose of this program is downloading the observation data of 160 stations 
# over Chine from China National Climate Center(NCC) 
#========================================================================

url=http://cmdp.ncc-cma.net/upload/upload8/t16
fig=/home/ys17-19/renql/project/2018homework/statistic/data/160preci/preci

#------begin downloading--------------------
for((i=1;i<=5;i++))
do
 #   echo "$fig"0"$i"
 #   wget -O "$fig"0"$i".txt "$url"0"$i"

# because the data is from 1951 Jan to 2018 Aug, so have delete the 2018 data from Jan to Aug 
    sed -i '537,544d' "$fig"0"$i".txt
done
 
#for((i=10;i<=12;i++))
#do
#    echo "$fig""$i"
#    wget -O "$fig""$i".txt "$url""$i"
#done


