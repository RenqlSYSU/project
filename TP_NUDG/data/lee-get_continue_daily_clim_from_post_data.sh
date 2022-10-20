#!/bin/sh
#-----------------------------------------------
#   This is a shell script for calculating the
# climatology from post-processed data, You 
# should set the basic parameters as below. 
# Good Luck!
#               Last Modified on  2015-09-21
#               Last Modified on  2017-04-03
#               A L_Zealot Product
#-----------------------------------------------

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
# With only ensemble output in it
#PRE_DIR_ORG=/users/yangsong3/L_Zealot/F/AMIP_C5PM/input_data2TP_NUDG_UVT
PRE_DIR_ORG=/users/yangsong3/renql/F/F2000_CAM5/ctrl_data/input2TP_NUDG

# Names of 2D fields
#FDNAME2D="(/\"FLUT\"/)" #often use
#FDNAME2D="(/\"TS\",\"PRECL\",\"PRECC\",\"PSL\",\"TMQ\"/)" #often use
FDNAME2D="(/\"PRECL\",\"PRECC\"/)" #often use

# Names of 3D fields
FDNAME3D="(/\"U\",\"V\",\"T\"/)" #often use
#FDNAME3D="(/\"V\",\"T\",\"OMEGA\",\"Q\",\"Z3\"/)" #often use
#FDNAME3D="(/\"U\",\"T\",\"Z3\"/)" #often use

# Names of 3D HY fields
#FDNAME3D_HY="(/\"U\",\"V\",\"T\",\"OMEGA\",\"Q\",\"RELHUM\",\"Z3\",\"DTCOND\"/)" # hybrid coordinate

# Process fig flag
FLAG2D=0
FLAG3D=1

#-----------------------------------------------------------
#Make dir
PRE_DIR=\"${PRE_DIR_ORG}\"
#mkdir -p $PRE_DIR_ORG/clim

#Output post processed 2D fields
if  [ $FLAG2D == 1 ] ; then
ncl -nQ pre_dir=$PRE_DIR            \
    fdname2d=$FDNAME2D          \
    ./extract-2D-clim-from-post-daily-esm.ncl
fi
#Output post processed 3D fields
if  [ $FLAG3D == 1 ] ; then
ncl -nQ pre_dir=$PRE_DIR            \
    fdname3d=$FDNAME3D          \
    ./extract-3D-clim-from-post-daily-esm.ncl
fi

