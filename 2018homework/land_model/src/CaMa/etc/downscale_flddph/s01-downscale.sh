#!/bin/sh

AREAS=`awk 'NR==2 {sub("area",""); print $0}' ./map/hires/location.txt`
echo $AREAS

for AREA in $AREAS
do

  FLDDPH="./out/flddph1990.mon"      # sample monthly flood depth file
  TREC=5                             # downscale May (irec=5)

  ./downscale_flddph $AREA $FLDDPH $TREC
done
