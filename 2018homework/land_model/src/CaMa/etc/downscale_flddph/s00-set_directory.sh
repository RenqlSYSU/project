#!/bin/sh

MAPDIR="../../map/global_15min/"
#OUTDIR="../../out/global_15min/"
OUTDIR="./sample/"                 ## global 15min sample flood depth file (monthly)

rm -f map
rm -f out
ln -s $MAPDIR map
ln -s $OUTDIR out

AREAS=`awk 'NR==2 {sub("area",""); print $0}' ./map/hires/location.txt`
echo $AREAS

