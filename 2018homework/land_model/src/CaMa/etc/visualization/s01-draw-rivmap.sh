#!/bin/bash

### set map dir ###
MAPDIR="../../map/global_15min/"
rm -f map
ln -s $MAPDIR map

### Plot File Name ##
PSFILE="rivmap.ps"

### Plot Domain ###
XMIN=-80                       # west
XMAX=-45                       # east
YMIN=-20                       # south
YMAX=5                         # north

### Plot Size ###
XSIZE=14                       # width of window
YSIZE=10                       # hight of window

AINT=5                         # Label Interbals
FINT=2.5                       # 

####################
GSIZE=`awk 'NR==5{print $1}' ./map/params.txt`

JFLAG="-JX${XSIZE}/${YSIZE}"
RFLAG="-R${XMIN}/${XMAX}/${YMIN}/${YMAX}"
BFLAG="-Ba${AINT}f${FINT}g0::/a${AINT}f${FINT}g0:::.:WneS"
FLAGS="-V -P ${JFLAG} ${RFLAG} ${BFLAG}"

RFLAG2="-R0/${XSIZE}/0/${YSIZE}"
BFLAG2="-Ba0f0g0:${XLABEL}:/a0f0g0:${YLABEL}::.${TITLE}:wnes"
FLAGS2="-V -P ${JFLAG} ${RFLAG2} ${BFLAG2}"
###################

psbasemap -K -X2 -Y2 $FLAGS -G255/255/255 > $PSFILE

### lsmask
./bin/print_grid $XMIN $XMAX $YMAX $YMIN $GNUM |\
xyz2grd -Gtmp.grd $RFLAG -I${GSIZE}/${GSIZE}
grdimage -O -K tmp.grd -Cgrid.cpt  $FLAGS >> $PSFILE

### flow direction
./bin/txt_vector $XMIN $XMAX $YMAX $YMIN > tmp.txt

for LEVEL in 01 02 03 04 05 06 07 08 09 10
do
  ./bin/print_rivvec tmp.txt 1 $LEVEL >  tmp2.txt
  awk '{print $1, $2, $4, $5}' tmp2.txt |\
  psxy -O -K ${FLAGS} -SvS0.${LEVEL}/0/0 -G50/50/255 >> $PSFILE
done

### bifurcation channel
#./bin/print_fldpth $XMIN $XMAX $YMAX $YMIN > fldpth.txt    ##  lon1, lat1, lon2,  lat2,  dst, elv, wth
#
#awk '$6>9000 {print $1, $2, $3, $4}' fldpth.txt |\
#psxy -O -K ${FLAGS} -SvS0.02/0/0 -G50/200/0 >> $PSFILE
#
#awk '$6<9000 {print $1, $2, $3, $4}' fldpth.txt |\
#psxy -O -K ${FLAGS} -SvS0.04/0/0 -G255/0/0 >> $PSFILE

#####
pstext -O ${FLAGS2} -N >> $PSFILE << EOF
0.5 10.5 20 0 1 5 River Network Map
EOF

#####

rm -f tmp.txt
rm -f tmp.grd
rm -f tmp2.txt
rm -f fldpth.txt
