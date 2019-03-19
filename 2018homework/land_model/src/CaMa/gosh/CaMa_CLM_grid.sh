#!/bin/sh
# =========================================================
# execute CaMa-Flood simulation
# =========================================================

CAS_RUNDIR=$1
CAS_OUTDIR=$2
CAS_RSTDIR=$3
DAT_RTMDIR=$4
TIMESTEP=$5
RUN_TYPE=$6
RESTART_FREQ=$7


# CaMa-Flood base directory
CAMADIR=`pwd`/..

cd ${CAMADIR}/gosh

# OPTIONS #
# You can change detailed setting by editing the namelist below.
# for entire list of options, see $CaMa-Flood/mod/mod_input.F

##### Basic Settings ##################
BASE=$CAMADIR                         # base directory

OUTDIR=$CAS_OUTDIR/CaMa               # CaMa output dir
mkdir -p $OUTDIR                      # make output dir
RESDIR=$CAS_RSTDIR/CaMa               # CaMa restart dir
mkdir -p $RESDIR                      # make restart dir

export OMP_NUM_THREADS=30             # OpenMP cpu num
LFLDOUT=".TRUE."                      # .TRUE. to activate floodplain discharge
LPTHOUT=".FALSE."                     # .TRUE. to activate bifurcation flow, mainly for delta simulation
LSTOONLY=".FALSE."                    # .TRUE. for restart only from storage (no previous-time-step discharge)


##### Time Step ####################### NOTE: spatial resolutions is set by dimension info file "diminfo.txt"
LADPSTP=".TRUE."                      # .TRUE. for adaptive time step
DT=$TIMESTEP                          # time step [sec]; set to 86400 for adaptive time step
DTIN=$TIMESTEP                        # input runoff time step [sec]


# !!NOTE: not applicable anymore, drived by CLM
##### Simulation Time #################
YSTART=1990                           # start year (from YSTART / Jan / 1st )
YEND=1991
SPINUP=2                              # 1 for restart, 2 for spinup
NSP=1                                 # spinup years

SPINUP=$RUN_TYPE                      # 1 for restart, 2 for spinup (initial run)

# CRESTSTO="set-by-shell"             # restart file name

##### Map & Topography ################
FMAP="${DAT_RTMDIR}"                  # map directory
#FMAP="${DAT_RTMDIR}/region_15min"    # (regional map) 

CDIMINFO="${FMAP}/diminfo_1deg.txt"   # dimention info (1deg, 0E->360E, 90N-90S)
#CDIMINFO="${FMAP}/diminfo_30min.txt" # dimention info (30min, 0E->360E, 90S->90N), generate a new matrix in map dir if needed

CNEXTXY=${FMAP}/nextxy.bin            # downstream xy (river network map)
CGRAREA=${FMAP}/grarea.bin            # unit-catchment area [m2]
CELEVTN=${FMAP}/elevtn.bin            # base elevation      [m]
CNXTDST=${FMAP}/nxtdst.bin            # downstream distance [m]
CRIVLEN=${FMAP}/rivlen.bin            # channel length      [m]
CFLDHGT=${FMAP}/fldhgt.bin            # floodplain elevation profile (height above 'elevtn') [m]

CRIVWTH=${FMAP}/rivwth.bin            # channel width       [m] (empirical power-low)
# CRIVWTH=${FMAP}/rivwth_gwdlr.bin    # channel width       [m] (GWD-LR + filled with empirical)
CRIVHGT=${FMAP}/rivhgt.bin            # channel depth       [m] (empirical power-low)
# CRIVHGT=${FMAP}/rivhgt_pth.bin      # channel depth       [m] (modified for bifurcation)

CPTHOUT=${FMAP}/fldpth.txt            # bifurcation channel list


##### Input Runoff Forcing ################# input runoff resolution should be consistent with "inpmat.bin"
LINTERP=".FALSE."                     # true: runoff interpolation using inpmat, false: nearest point interpolation
CINPMAT=${FMAP}/inpmat-1deg.bin       # runoff input matrix (1deg, 0E->360E, 90N-90S)     !! for sample bonary input

#--------------------------------------
# CINPMAT=${FMAP}/inpmat-30min.bin    # runoff input matrix (30min, 0E->360E, 90S->90N) !! for sample netCDF input
#                                     # generate a new matrix in map dir if needed
#LINPCDF=".FALSE."                    # true: netCDF input file
#CROFDIR="${BASE}/inp/ELSE_GPCC/Roff/"# runoff directory
#CROFPRE="Roff____"                   # runoff prefix/suffix  ( $(PREFIX)yyyymmdd$(SUFFIX) )
#CROFSUF=".one"
#
#LINPCDF=".TRUE."
#CROFDIR="${BASE}/inp/ELSE_GPCC/runoff_nc/" #   runoff directory
#CROFCDF="set-by-shell"               # netCDF runoff filename
#CROFVAR="runoff"                     # netCDF runoff variable name
#SYEARIN="set-by-shell"               # netCDF runoff file, start date
#SMONIN="set-by-shell"
#SDAYIN="set-by-shell"
#--------------------------------------

DROFUNIT=1.D-3                        # runoff unit conversion (1.D-3 when input [mm] is converted to [m3/m2])

##### Output Settings #################
LOUTCDF=".FALSE."                     # true for netCDF output, false for plain binary output
COUTDIR=""                            # output directory

# output variables set "NONE" for no output
CRIVOUTDIR="NONE"                     # river discharge         [m3/s]
CRIVSTODIR="NONE"                     # river storage           [m3]
CRIVVELDIR="NONE"                     # river flow velocity     [m/s]
CRIVDPHDIR="${OUTDIR}/"                       # river water depth       [m]

CFLDOUTDIR="NONE"                     # floodplain discharge    [m3/s]
CFLDSTODIR="NONE"                     # floodplain storage      [m]
CFLDDPHDIR="${OUTDIR}/"               # floodplain water depth  [m]
CFLDAREDIR="${OUTDIR}/"               # flooded area            [m]
CFLDFRCDIR="${OUTDIR}/"               # flooded area fraction   [m2/m2]

CSFCELVDIR="${OUTDIR}/"               # water surface elevation [m]
COUTFLWDIR="${OUTDIR}/"               # total discharge (rivout+fldout)   [m3/s]
CSTORGEDIR="${OUTDIR}/"               # total storage (rivsto+fldsto)     [m3]

CPTHOUTDIR="${OUTDIR}/"               # net bifurcation flow (grid-based) [m3/s]
CPTHFLWDIR="${OUTDIR}/"               # bifurcation flow (channel based)  [m3/s]

##### Model Parameters ################
PMANRIV=0.03D0                        # manning coefficient river
PMANFLD=0.10D0                        # manning coefficient floodplain
PCADP=0.7                             # satety coefficient for CFL condition

##### Spatial Resolutions #############
# Set by "diminfo.txt".
# For manual setting, set CDIMINFO="NONE"
# NX="set-by-diminfo"                 # number of grids in east-west
# NX="set-by-diminfo"                 # number of grids in east-west
# NLFP="set-by-diminfo"               # floodplain layer

# NXIN="set-by-diminfo"               # number of input grids in east-west
# NYIN="set-by-diminfo"               # number of input grids in east-west
# INPN="set-by-diminfo"               # max number of input grids for one cama grid

# WEST="set-by-diminfo"               # domain west  edge
# EAST="set-by-diminfo"               # domain east  edge
# NORTH="set-by-diminfo"              # domain north edge
# SOUTH="set-by-diminfo"              # domain south edge

### End of Setting ####################


## create running dir 
cd ${CAS_RUNDIR}

## if new simulation, remove old files in running directory

if [ $SPINUP -eq 2 ]; then
  rm -f ${CAS_RUNDIR}/????-sp*
  rm -f ${CAS_RUNDIR}/*.bin
  rm -f ${CAS_RUNDIR}/*.pth
  rm -f ${CAS_RUNDIR}/*.vec
  rm -f ${CAS_RUNDIR}/*.nc
  rm -f ${CAS_RUNDIR}/*.txt
  rm -f ${CAS_RUNDIR}/restart*.*
fi

if [ $SPINUP -eq 2 ];then
  IRESTART=2
  CRESTSTO=""
else
  IRESTART=1
  CRESTSTO=""
fi

# set restart dir
########### create input namelist ##########
cat > ${CAS_RUNDIR}/input_flood.nam << EOF
&NRUNVER
IRESTART=$IRESTART                  ! 1=> restart;  2=>spinup
CRESTDIR="$RESDIR/"                 ! restart directory
CRESTSTO="$CRESTSTO"                ! restart file
LSTOONLY=$LSTOONLY                  ! true for restart only from storage
LRESTCDF=.FALSE.                    ! true for netCDF restart file
RESTFREQ=$RESTART_FREQ              ! 0: yearly restart file, 1: daily restart file, 2: monthly restart file
/
&NSIMTIME
ISYEAR=$YSTART                      ! start year
ISMON=1                             ! month 
ISDAY=1                             ! day        (assumed at 00UTC)
IEYEAR=$YEND                        ! end year
IEMON=1                             ! end
IEDAY=1                             ! end        (assumed at 00UTC)
/
&NMAP
LMAPCDF=.false.                     ! true for netCDF map input
CDIMINFO="${CDIMINFO}"              ! dimention info
CNEXTXY="${CNEXTXY}"                ! downstream xy (river network map)
CGRAREA="${CGRAREA}"                ! unit-catchment area [m2]
CELEVTN="${CELEVTN}"                ! base elevation      [m]
CNXTDST="${CNXTDST}"                ! downstream distance [m]
CRIVWTH="${CRIVWTH}"                ! channel width       [m]
CRIVLEN="${CRIVLEN}"                ! channel length      [m]
CRIVHGT="${CRIVHGT}"                ! channel depth       [m]
CFLDHGT="${CFLDHGT}"                ! floodplain elevation profile [m]
CPTHOUT="${CPTHOUT}"                ! bifurcation channel list
CRIVCLINC="NONE"                    ! * netCDF river maps
CRIVPARNC="NONE"                    ! * netCDF river width & depth
/
&NINPUT 
LINTERP=${LINTERP}                  ! true for runoff interpolation using input matrix
LINPCDF=${LINPCDF}                  ! true for netCDF input
CINPMAT="${CINPMAT}"                ! input matrix file name
CRUNOFFDIR="${CROFDIR}"             ! runoff input directory
CRUNOFFPRE="${CROFPRE}"             ! runoff input prefix
CRUNOFFSUF="${CROFSUF}"             ! runoff input suffix
CRUNOFFCDF="${CROFCDF}"             ! * netCDF input runoff file name
CROFCDFVAR="${CROFVAR}"             ! * netCDF input runoff variable name
SYEARIN=$SYEARIN                    ! * for netCDF input start date (start of the initial time step)
SMONIN=12
SDAYIN=31
LINTERPCDF=.FALSE.                  ! * true for netCDF input matrix
/
&NOUTPUT
LOUTCDF=${LOUTCDF}                  ! true for netCDF output
COUTDIR="${COUTDIR}"                ! output directory ("NONE" for no output)
CRIVOUTDIR="${CRIVOUTDIR}"          ! river discharge        [m3/s]
CRIVSTODIR="${CRIVSTODIR}"          ! river storage          [m3]
CRIVVELDIR="${CRIVVELDIR}"          ! river flow velocity    [m/s]
CRIVDPHDIR="${CRIVDPHDIR}"          ! river water depth      [m]
CFLDOUTDIR="${CFLDOUTDIR}"          ! floodplain discharge   [m3/s]
CFLDSTODIR="${CFLDSTODIR}"          ! floodplain storage     [m3]
CFLDDPHDIR="${CFLDDPHDIR}"          ! floodplain water depth [m]
CFLDFRCDIR="${CFLDFRCDIR}"          ! flooded area fraction  [m2/m2]
CFLDAREDIR="${CFLDAREDIR}"          ! flooded area           [m2]
CSFCELVDIR="${CSFCELVDIR}"          ! water surface elevation           [m]
COUTFLWDIR="${COUTFLWDIR}"          ! total discharge (rivout+fldout)   [m3/s]
CSTORGEDIR="${CSTORGEDIR}"          ! total storage   (rivsto+fldsto)   [m3]
CPTHOUTDIR="${CPTHOUTDIR}"          ! net bifurcation flow (grid-based) [m3/s]
CPTHFLWDIR="${CPTHFLWDIR}"          ! bifurcation flow (channel-based)  [m3/s]
COUTINSDIR="NONE"                   ! instantaneous discharge (no river routing, summation of upstream runoff)
LOUTCDF=.TRUE.                      ! output in NetCDF format 
LOUTVEC=.FALSE.                     ! for 1-D land-only output (small data size, post processing required)
/
&NCONF                              ! * NX, NY, NFLP, NXIN, NYIN, INPN, WEST, EAST, NORTH, SOUTH set by diminfo.txt
DT=$DT                              ! time step [sec]
DTIN=$DTIN                          ! input runoff time step [sec]
DROFUNIT=$DROFUNIT                  ! runoff unit conversion (1.D-3 when input [mm] is converted to [m3/m2]
LADPSTP=$LADPSTP                    ! true for adaptive time step
LFLDOUT=$LFLDOUT                    ! true to activate floodplain discharge
LPTHOUT=$LPTHOUT                    ! true to activate bifurcation channel flow
LFLD=.TRUE.                         ! true to activate floodplain inundation
LKINE=.FALSE.                       ! true for kinematic river routing
LMAPEND=.FALSE.                     ! true to convert map data endian
LINPEND=.FALSE.                     ! true to convert input data endian
LLEAPYR=.TRUE.                      ! true for leap year calculatuon, false: always 365days/year
/
&NPARAM
PMANRIV=$PMANRIV                    ! manning coefficient river
PMANFLD=$PMANFLD                    ! manning coefficient floodplain
PGRV=9.8D0                          ! accerelation due to gravity
PDSTMTH=10000.D0                    ! downstream distance at river mouth [m]
PCADP=$PCADP                        ! satety coefficient for CFL condition
PMINSLP=1.D-5                       ! * minimum slope (for kinematic wave)
/
EOF

exit 0
