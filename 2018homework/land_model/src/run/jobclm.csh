#!/bin/csh

#-------------------------------------------------------
# [1] Define the JOB
#-------------------------------------------------------

set RUN_CLM_SRF="NO"     	# "YES" = MAKE CoLM surface characteristic data
                                # "NO"  = NOT make CoLM surface characteristic data

set RUN_CLM_INI="YES"    	# "YES' = MAKE CoLM initial data
                                # "NO"  = Restart run

set RUN_CaMa="NO"       	# "YES" = OPEN CaMa-Flood
                                # "NO"  = CLOSE CaMa-Flood [No river routing]

set RUN_CLM="YES"        	# "YES" = RUN CoLM
                                # "NO'  = NOT RUN CoLM


# case name and simulating time setting
#-------------------------------------------------------
set CASE_NAME   = CLMCRU           	# case name                                            <MARK #1>
set GREENWICH   = .true.        	# 'true' for greenwich time, 'false' for local time
set START_YEAR  = 2001          	# model start year                                     <MARK #2>
set START_DAY   = 1             	# model start Julian day
set START_SEC   = 0             	# model start secs of day
set END_YEAR    = 2003          	# model end year
set END_DAY     = 1             	# model end Julian day
set END_SEC     = 0             	# model end secs of day
set SPINUP_YEAR = $START_YEAR   	# spin-up end year, set default to SATRT_YEAR
set SPINUP_DAY  = $START_DAY    	# spin-up end day, set default to START_DAY
set SPINUP_SEC  = $START_SEC    	# spin-up end sec, set default to START_SEC
set TIMESTEP    = 1800.         	# model time step

set WOUT_FREQ   = MONTHLY         	# write output  file frequency: HOURLY/DAILY/MONTHLY/YEARLY
set WRST_FREQ   = MONTHLY		# write restart file frequency: HOURLY/DAILY/MONTHLY/YEARLY

# model resolution and running scope setting
#-------------------------------------------------------
set LON_POINTS  =  720              
set LAT_POINTS  =  360               
set EDGE_N      =   90.
set EDGE_E      =  180.
set EDGE_S      =  -90.
set EDGE_W      = -180.

# set forcing observational height (unit: m)
#-------------------------------------------------------
set HEIGHT_V    = 100.
set HEIGHT_T    =  50.
set HEIGHT_Q    =  50.


#-------------------------------------------------------
# [2] Set necessary environment variables
#-------------------------------------------------------

# clm src directory
#-------------------------------------------------------
setenv CLM_ROOT   $HOME/clm                                          # <MARK #3>
setenv CLM_INCDIR $CLM_ROOT/include
setenv CLM_SRFDIR $CLM_ROOT/mksrfdata
setenv CLM_INIDIR $CLM_ROOT/mkinidata
setenv CLM_SRCDIR $CLM_ROOT/main

# inputdata directory
setenv DAT_ROOT   /p13/data/inputdata                                # <MARK #4>
setenv DAT_RAWDIR /p13/data/CLMrawdata         
setenv DAT_ATMDIR $DAT_ROOT/atm/cruncep
setenv DAT_SRFDIR $DAT_ROOT/srf/global_0.5x0.5
setenv DAT_RTMDIR $DAT_ROOT/rtm/global_15min

# case directory
#-------------------------------------------------------
setenv CAS_ROOT   /p13/cases                                         # <MARK #5>
setenv CAS_RUNDIR $CAS_ROOT/$CASE_NAME
setenv CAS_OUTDIR $CAS_RUNDIR/output
setenv CAS_RSTDIR $CAS_RUNDIR/restart

mkdir -p $DAT_SRFDIR
mkdir -p $CAS_RUNDIR
mkdir -p $CAS_OUTDIR
mkdir -p $CAS_RSTDIR




#------------------------------------------------------
# [3] build define.h in ./include directory
#------------------------------------------------------

\cat >! .tmp << EOF
#define	USGS_CLASSIFICATION       ! USGS/IGBP
#undef	RDGRID                    !
#undef	RAWdata_update            !
#undef	DYN_PHENOLOGY             !
#undef	SOILINI                   !
#define	LANDONLY                  !
#undef	LAND_SEA                  !
#undef	SOIL_REFL_GUESSED         !
#define	SOIL_REFL_READ            !
#define	WO_${WOUT_FREQ}           !
#define	WR_${WRST_FREQ}           !
#undef	CLMDEBUG                  !
#define	USE_CRUNCEP_DATA          ! QIAN/PRINCETON/CRUNCEP/POINT
#define	HEIGHT_V $HEIGHT_V        ! 
#define	HEIGHT_T $HEIGHT_T        !
#define	HEIGHT_Q $HEIGHT_Q        !
#define	OPENMP 28                 !
EOF

#-------------------------------------------------------#
#              --- USER SETTING END ---                 #
# DO NOT EDIT THE BELOW SCRIPTS UNLESS YOU KNOW EXACTLY #
# WHAT YOU ARE DOING                                    #
#-------------------------------------------------------#

sed -i 's/\!.*//g' .tmp
sed -i '/^ *$/d' .tmp
mv -f .tmp $CLM_INCDIR/define.h

if ( $RUN_CaMa == "YES" ) then
  echo "#define CaMa_Flood" >> $CLM_INCDIR/define.h
else
  echo "#undef  CaMa_Flood" >> $CLM_INCDIR/define.h
endif



#-------------------------------------------------------
# [4] compling and executing CoLM surface data making
#-------------------------------------------------------
if ( $RUN_CLM_SRF == "YES" ) then

# Compile
cd $CLM_SRFDIR
make clean
make >& $CAS_RUNDIR/compile.mksrf.log || exit 5

# Create an input parameter namelist file
\cat >! $CAS_RUNDIR/mksrf.stdin << EOF
&mksrfexp
dir_rawdata        = '$DAT_RAWDIR/'
dir_model_landdata = '$DAT_SRFDIR/'
lon_points         = $LON_POINTS
lat_points         = $LAT_POINTS
edgen              = $EDGE_N
edgee              = $EDGE_E
edges              = $EDGE_S
edgew              = $EDGE_W
/
EOF

# Executing CoLM initialization'
$CLM_SRFDIR/srf.x < $CAS_RUNDIR/mksrf.stdin > $CAS_RUNDIR/exe.mksrf.log || exit 4

echo 'Making the CoLM surface data completed'

endif



#-------------------------------------------------------
# [5] compling and executing CoLM initialization
#-------------------------------------------------------

cd $CLM_INIDIR

# Create an input parameter namelist file
\cat >! $CAS_RUNDIR/mkini.stdin << EOF
&clminiexp
site = '$CASE_NAME'
dir_model_landdata = '$DAT_SRFDIR/'
dir_restart_hist   = '$CAS_RSTDIR/'
dir_infolist       = '$CAS_RUNDIR/'
lon_points         = $LON_POINTS
lat_points         = $LAT_POINTS
greenwich          = $GREENWICH
s_year             = $START_YEAR
s_julian           = $START_DAY
s_seconds          = $START_SEC
/
EOF

if ( $RUN_CLM_INI == "YES" ) then

# CoLM initialization for startup run
#-------------------------------------------------------
make clean
make >& $CAS_RUNDIR/compile.mkini.log || exit 5

$CLM_INIDIR/initial.x < $CAS_RUNDIR/mkini.stdin > $CAS_RUNDIR/exe.mkini.log || exit 4
echo 'CoLM initialization completed'

else

# for restart run
#-------------------------------------------------------
if (! -e $CAS_RUNDIR/clmini.infolist) then
  echo 'ERROR: no initial run detected, please run clm initialization first!'; exit
endif

sed -e    "s/s_year *=.*/s_year    = ${START_YEAR}/"  \
    -e  "s/s_julian *=.*/s_julian  = ${START_DAY}/"   \
    -e "s/s_seconds *=.*/s_seconds = ${START_SEC}/"   \
< $CAS_RUNDIR/clmini.infolist > .tmp

mv -f .tmp $CAS_RUNDIR/clmini.infolist

echo 'CoLM initialization for restart run completed'

endif



#-------------------------------------------------------
# [6] compliling CaMa Flood Model and make namelist file
#-------------------------------------------------------
if ( $RUN_CaMa == "YES" ) then

echo 'Compiling and initilizing CaMa'

setenv CaMa_DIR $CLM_ROOT/CaMa

set RESTART = 1
set SPINUP  = 2

set RESTART_FREQ = 2
if ( $WRST_FREQ == "YEARLY"  ) then 
  set RESTART_FREQ = 0 
endif
if ( $WRST_FREQ == "DAILY"   ) then 
  set RESTART_FREQ = 1 
endif
if ( $WRST_FREQ == "MONTHLY" ) then 
  set RESTART_FREQ = 2 
endif

# compile
cd $CaMa_DIR/gosh
chmod u+x compile.sh
./compile.sh >& $CAS_RUNDIR/compile.CaMa.log || exit 5
echo 'Compiling CaMa Flood Model completed'

# Create an input parameter namelist file for CaMa Flood Model
chmod u+x CaMa_CLM_grid.sh

if ($RUN_CLM_INI == "YES") then
./CaMa_CLM_grid.sh $CAS_RUNDIR $CAS_OUTDIR $CAS_RSTDIR $DAT_RTMDIR $TIMESTEP $SPINUP  $RESTART_FREQ
else
./CaMa_CLM_grid.sh $CAS_RUNDIR $CAS_OUTDIR $CAS_RSTDIR $DAT_RTMDIR $TIMESTEP $RESTART $RESTART_FREQ
endif

echo 'CaMa compiling and initialization completed'

endif



#-------------------------------------------------------
# [7] compiling and executing CoLM model
#-------------------------------------------------------
if ( $RUN_CLM == "YES" ) then

# Compile
cd $CLM_SRCDIR
make clean
rm -f $CAS_RUNDIR/compile.main.log

make >& $CAS_RUNDIR/compile.main.log || exit 5

echo 'Compiling CoLM completed'

# Create an input parameter namelist file
\cat >! $CAS_RUNDIR/timeloop.stdin << EOF
&clmexp
site = '$CASE_NAME'
dir_model_landdata = '$DAT_SRFDIR/'
dir_forcing        = '$DAT_ATMDIR/'
dir_output         = '$CAS_OUTDIR/'
dir_restart_hist   = '$CAS_RSTDIR/'
lon_points         = $LON_POINTS
lat_points         = $LAT_POINTS
deltim             = $TIMESTEP
solarin_all_band   = .true.
e_year             = $END_YEAR
e_julian           = $END_DAY
e_seconds          = $END_SEC
p_year             = $SPINUP_YEAR
p_julian           = $SPINUP_DAY
p_seconds          = $SPINUP_SEC
EOF

\cat $CAS_RUNDIR/clmini.infolist >> $CAS_RUNDIR/timeloop.stdin

# Executing the CoLM
echo 'Executing CoLM...'

#----------------------------------------------------------------------
cd $CAS_RUNDIR

rm -f $CAS_RUNDIR/exe.timeloop.log
ln -sf $CLM_SRCDIR/clm.x $CAS_RUNDIR/.

/usr/bin/time ./clm.x < $CAS_RUNDIR/timeloop.stdin > $CAS_RUNDIR/exe.timeloop.log || exit 4

#----------------------------------------------------------------------

endif

echo 'CoLM Execution Completed'


echo '-----------------------------------------------------------------'
echo ' End of c-shell script                                           '
echo '-----------------------------------------------------------------'
