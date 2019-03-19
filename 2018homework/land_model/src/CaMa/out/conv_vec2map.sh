#!/bin/sh
# convert 1D land only output to 2D map

OUTDIR=$1
OUTDIR=`echo $OUTDIR | tr -d "\/"`
#OUTDIR="global_15min"

SYEAR=1990
EYEAR=1991

IYEAR=$SYEAR
while [ $IYEAR -le $EYEAR ];
do

#  VARS='rivout rivsto rivvel rivdph fldsto flddph fldfrc fldare sfcelv'
  VARS="rivout"
  for VAR in $VARS
  do
    echo "./vec2map ${OUTDIR}/ $VAR $IYEAR"
    ./vec2map ${OUTDIR}/ $VAR $IYEAR
  done

IYEAR=`expr $IYEAR + 1`
done

