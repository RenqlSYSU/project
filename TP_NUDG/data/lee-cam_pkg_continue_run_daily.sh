#!/bin/bash
#-----------------------------------------------
#   This is a shell script for configuring the basic post processing tool of CAM model, You 
# should set the basic parameters as below. Good Luck!
#               Last Modified on  2016-04-02
#               A L_Zealot Product
#Modified by Ql_Ren 2017/12/26 for TP-NUDG
#-----------------------------------------------

# Path of the original data
# Caution: DO NOT DELETE \" IN STRING!
PRE_DIR[1]=\"/home/ys17-23/Extension2/renql/model_data/F2000_CAM5/pre_ctrl/\"

# Path of the post processed data
PRO_DIR[1]=\"/home/ys17-23/Extension2/renql/model_data/F2000_CAM5/\"

# Case name
CASENAME[1]=\"F2000_CAM5\"

# Names of 2D fields
#FDNAME2D="(/\"PRECL\",\"PRECC\"/)" #often use
#FDNAME2D="(/\"PRECL\",\"PRECC\",\"LHFLX\",\"PS\",\"PSL\",\"QFLX\",\"TS\",\"TMQ\"/)" #often use
FDNAME2D="(/\"PRECL\",\"PRECC\",\"PS\",\"PSL\",\"TS\"/)" #often use
#FDNAME2D="(/\"PSL\"/)" #often use
#FDNAME2D="(/\"PRECT\",\"U850\",\"V850\"/)" #often use
#FDNAME2D="(/\"LHFLX\",\"SHFLX\"/)" #often use

# Names of 3D fields
#FDNAME3D="(/\"U\",\"V\"/)" # hybrid coordinate
#FDNAME3D="(/\"U\",\"V\",\"T\"/)" #often use
FDNAME3D="(/\"U\",\"V\",\"T\",\"OMEGA\",\"Z3\",\"Q\"/)" # hybrid coordinate
#FDNAME3D="(/\"Q\"/)" # hybrid coordinate
#FDNAME3D_HY="(/\"U\",\"V\",\"T\",\"OMEGA\",\"Q\",\"RELHUM\",\"Z3\",\"DTCOND\"/)" # hybrid coordinate
#FDNAME3D_HY="(/\"U\",\"V\",\"T\"/)" # hybrid coordinate

# First year of the subset
FRSTYEAR=0001 #1979 #

# Last year of the subset
LSTYEAR=0030 #2005 #

# Layers of 3D fields
# CAM4 = 26; CAM5 = 30
LAYERS=30

# Output specific pressure layers
# CAUTION: Do not leave species between element!
#PLEV="(/925,850,700,600,500,400,300,200,100,50/)" #10levels
#PLEV="(/1000,925,850,700,500,200/)"
#PLEV="(/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/)" #18levels
#PLEV="(/1000,925,850,700,600,500,400,300/)" #used for Q
PLEV="(/1000,925,850,700,600,500,400,350,300,250,200,150,100/)" #17levels,50,20,10,5
#PLEV="(/1000,925,850,700,600,500,400,350,300,250,200/)" #14levels

extrap_option=1
# 1 mean Use the vinth2p_ecmwf for the extrapolation below ground.
# 0 mean no extrapolation when the pressure level is outside of the range of psfc.
# A scalar integer indicating which variable to interpolate for vinth2p_ecmwf
# 1 = temperature, -1 = geopotential height, 0 = all others.
varflg="(/0,0,1,0,-1,0/)"

# Process flag
FLAG_2D=0
FLAG_3D=1
FLAG_3D_HY=0

#-----------------------------------------------------------

for ni in {1..1}
do
echo "${CASENAME[ni]}"
echo "${PRO_DIR[ni]}"
#Output post processed 2D fields
if [ $FLAG_2D == 1 ] ; then
    echo "-----package 2D field    (1)-----"
    ncl -nQ \
        pre_dir=${PRE_DIR[ni]}            \
        pro_dir=${PRO_DIR[ni]}            \
        fdname2d=$FDNAME2D          \
        frstyear=$FRSTYEAR          \
        lstyear=$LSTYEAR           \
        case_name=${CASENAME[ni]}         \
        ./package_2D_from_raw_data_daily-160402.ncl
fi


#Output post processed 3D fields
if [ $FLAG_3D == 1 ] ; then
    echo "-----package 3D field    (1)-----"
    ncl -nQ \
        pre_dir=${PRE_DIR[ni]}            \
        pro_dir=${PRO_DIR[ni]}            \
       fdname3d=$FDNAME3D          \
       layers=$LAYERS              \
       plev=$PLEV                  \
       frstyear=$FRSTYEAR          \
       lstyear=$LSTYEAR           \
        case_name=${CASENAME[ni]}         \
        extrap_option=$extrap_option         \
        varflg=$varflg \
       ./package_3D_from_raw_data_daily-160808.ncl
fi

#Output post processed 3D Hybrid fields
if [ $FLAG_3D_HY == 1 ] ; then
    echo "-----package 3D field    (1)-----"
    ncl -nQ \
        pre_dir=${PRE_DIR[ni]}            \
        pro_dir=${PRO_DIR[ni]}            \
       fdname3d=$FDNAME3D_HY          \
       frstyear=$FRSTYEAR          \
       lstyear=$LSTYEAR           \
       layers=$LAYERS              \
        case_name=${CASENAME[ni]}         \
       ./package_3DHY_from_raw_data_daily-171019.ncl
fi
done
