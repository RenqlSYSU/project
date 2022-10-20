#!/bin/sh
#-----------------------------------------------
#   This is a shell script for configuring the basic post processing tool of CAM model, You 
# should set the basic parameters as below. Good Luck!
#               Last Modified on  2016-04-02
#               A L_Zealot Product
#Modified by Ql_Ren 2017/12/26 for TP-NUDG
#-----------------------------------------------

# Path of the original data
# Caution: DO NOT DELETE \" IN STRING!
#PRE_DIR[1]=\"/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/\"
#PRE_DIR[2]=\"/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/pre_data_24h/\"
#PRE_DIR[3]=\"/users/yangsong3/L_Zealot/F/AMIP_C5PM/exe/\"
#PRE_DIR[1]=\"/home/ys17-19/renql/model/TP_CR/\"
#PRE_DIR[2]=\"/home/ys17-19/renql/model/TP_CTRL/pre/\"
PRE_DIR[1]=\"/users/yangsong3/renql/F/F2000_CAM5/ctrl_data/\"

# Path of the post processed data
#PRO_DIR[1]=\"/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/post_data_6h/\"
#PRO_DIR[2]=\"/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/post_data_24h/\"
#PRO_DIR[3]=\"/users/yangsong3/L_Zealot/F/AMIP_C5PM/post_data/\"
#PRO_DIR[1]=\"/home/ys17-19/renql/model/TP_CR/pro/\"
#PRO_DIR[2]=\"/home/ys17-19/renql/model/TP_CTRL/pro/\"
PRO_DIR[1]=\"/users/yangsong3/renql/F/F2000_CAM5/ctrl_data/input2TP_NUDG/\"

# Case name
#CASENAME[1]=\"AMIP_C5PM_TP_NUDG\"
#CASENAME[2]=\"AMIP_C5PM_TP_NUDG\"
#CASENAME[3]=\"AMIP_C5PM\"
#CASENAME[1]=\"TP_CR\"
#CASENAME[2]=\"TP_CTRL\"
CASENAME[1]=\"F2000_CAM5\"

# Names of 2D fields
#FDNAME2D="(/\"PRECL\",\"PRECC\",\"LHFLX\",\"PS\",\"PSL\",\"QFLX\",\"TS\",\"TMQ\"/)" #often use
FDNAME2D="(/\"PS\"/)" #often use
#FDNAME2D="(/\"PRECT\",\"U850\",\"V850\"/)" #often use
#FDNAME2D="(/\"PRECC\",\"PRECL\"/)" #often use
#FDNAME2D="(/\"LHFLX\",\"SHFLX\"/)" #often use

# Names of 3D fields
FDNAME3D="(/\"U\",\"V\"/)" # hybrid coordinate
#FDNAME3D="(/\"U\",\"V\",\"T\",\"OMEGA\",\"Q\",\"Z3\"/)" #often use
#FDNAME3D="(/\"DTCOND\"/)" #often use
#FDNAME3D_HY="(/\"RELHUM\"/)" #often use
#FDNAME3D_HY="(/\"U\",\"V\",\"T\",\"OMEGA\",\"Q\",\"RELHUM\",\"Z3\",\"DTCOND\"/)" # hybrid coordinate
FDNAME3D_HY="(/\"U\",\"V\",\"T\"/)" # hybrid coordinate

# First year of the subset
FRSTYEAR=0001

# Last year of the subset
LSTYEAR=0030

# Layers of 3D fields
# CAM4 = 26; CAM5 = 30
LAYERS=30

# Output specific pressure layers
# CAUTION: Do not leave species between element!
#PLEV="(/925,850,700,600,500,400,300,200,100,50/)"
#PLEV="(/1000,925,850,700,500,200/)"
PLEV="(/1000,925,850,700,600,500,400,300,200/)"

# Process flag
FLAG_2D=0
FLAG_3D=0
FLAG_3D_HY=1

#-----------------------------------------------------------

for ni in {1..1}
do
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
