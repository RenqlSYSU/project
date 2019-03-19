#!/bin/sh
# convert daily output to monthly average

OUTDIR=$1
OUTDIR=`echo $OUTDIR | tr -d "\/"`
#OUTDIR="global_15min"

SYEAR=1990
EYEAR=1990

IYEAR=$SYEAR
while [ $IYEAR -le $EYEAR ];
do
#  VARS='rivout rivsto rivvel rivdph fldout fldsto flddph fldfrc fldare sfcelv outflw storge'
  VARS='flddph'
  for VAR in $VARS
  do
    echo "./day2mon ${OUTDIR}/ $VAR $IYEAR"
    ./day2mon ${OUTDIR}/ $VAR $IYEAR
  done

IYEAR=`expr $IYEAR + 1`
done

