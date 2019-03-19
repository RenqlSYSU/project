module init_inputnam_mod
contains 
subroutine init_inputnam
!==================================================
! set up of namelist inputs.
!==================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: nsetfile, nulnam,  tmpnam, &
                      nrunver,  irestart,crestdir,creststo,lrestcdf,lstoonly,restfreq,&
                      nsimtime, isyear,  ismon,   isday,   ieyear,  iemon,   ieday,&
                      nmap,     cdiminfo,cnextxy, cgrarea, celevtn, cnxtdst, crivwth, crivlen, crivhgt, cfldhgt,&
                                cpthout, crivclinc,crivparnc,lmapcdf, &
                      ninput,   linterp, linpcdf, cinpmat,   crunoffdir,crunoffpre,crunoffsuf,crunoffcdf,crofcdfvar, &
                      noutput,  coutdir, crivoutdir,crivstodir,crivveldir,crivdphdir,&
                                         cfldoutdir,cfldstodir,cflddphdir,cfldfrcdir,cfldaredir,&
                                         csfcelvdir,coutflwdir,cstorgedir,cpthoutdir,cpthflwdir,&
                                         coutinsdir,loutcdf,   loutvec,   loutyyyy, &
                      nconf,    nx,ny,nlfp,dt,dt_def, nxin,nyin,dtin,drofunit,inpn,west,east,north,south,&
                                syearin, smonin,  sdayin, lfld,lkine,ladpstp,lfldout, lpthout, ldam, &
                                lmapend, linpend, lleapyr, &
                      nparam,   pmanriv, pmanfld, pgrv,   pdstmth,   pcadp,  pminslp, linterpcdf

use mod_map    ,only: regionthis
use mod_output ,only: coutvardir

implicit none

! ==============================
! *** 0. set input unit and open file 
open(nsetfile,file="input_flood.nam",status="old")

if(regionthis==1) then
  write(nulnam,*) " "
  write(nulnam,*) "init_inputnam: input_flood.nam open in unit:",nsetfile 
endif
! ==============================
! *** 1. basic simulation setting
! * defaults
irestart=1                       !! 1: start from "zero-storage", 2: restart from "restart file"
crestdir="./"                    !! restart file directory
creststo="restartin"             !! restart file name
lrestcdf=.false.                 !! true for netcdf restart file
lstoonly=.false.                 !! true for restart only from storage (no discharge restart, not strict)
restfreq=0                       !! 0: restart file at end of simulation, 1: daily restart file

! * change
read(nsetfile,nrunver)
if(regionthis==1) then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 1. init_runver****"
  write(nulnam,*) "irestart=",irestart
  write(nulnam,*) "crestdir=",trim(crestdir)
  write(nulnam,*) "lrestcdf=",lrestcdf
  write(nulnam,*) "lstoonly=",lstoonly
  write(nulnam,*) "restfreq=",restfreq
endif


! ==============================
! *** 2. set simulation time 
! * defaults
isyear=1990
ismon=1
isday=1
ieyear=1991
iemon=1
ieday=1

! * change
rewind(nsetfile)
read(nsetfile,nsimtime)
if(regionthis==1) then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 2. set_simtime****"
  write(nulnam,*) "isyear,ismon,isday:",isyear,ismon,isday
  write(nulnam,*) "ieyear,iemon,ieday:",ieyear,iemon,ieday
endif

! ==============================
! *** 3. set map
! * defaults
cdiminfo="NONE"
cnextxy="./nextxy.bin "
cgrarea="./grarea.bin"
celevtn="./elevtn.bin"
cnxtdst="./nxtdst.bin "
crivwth="./rivwth.bin"
crivlen="./rivlen.bin"
crivhgt="./rivhgt.bin"
cfldhgt="./fldhgt.bin"
cpthout="./fldpath.txt"
crivclinc="NONE"
crivparnc="NONE"
lmapcdf=.false.

! * change
rewind(nsetfile)
read(nsetfile,nmap)
if(regionthis==1) then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 3. set_setmap****"
  write(nulnam,*) "cdiminfo:",trim(cdiminfo)
  write(nulnam,*) "cnextxy:",trim(cnextxy)
  write(nulnam,*) "cgrarea:",trim(cgrarea)
  write(nulnam,*) "celevtn:",trim(celevtn)
  write(nulnam,*) "cnxtdst:",trim(cnxtdst)
  write(nulnam,*) "crivwth:",trim(crivwth)
  write(nulnam,*) "crivlen:",trim(crivlen)
  write(nulnam,*) "crivhgt:",trim(crivhgt)
  write(nulnam,*) "cfldhgt:",trim(cfldhgt)
  write(nulnam,*) "cpthout:",trim(cpthout)

  write(nulnam,*) "lmapcdf:"  ,lmapcdf
  write(nulnam,*) "crivclinc:",trim(crivclinc)
  write(nulnam,*) "crivparnc:",trim(crivparnc)
endif

! ==============================
! *** 4. set input
! * defaults
linterp=.false.                 !! use nearest interpolation (.true. for interpolation with input matrix)
cinpmat="NONE"

linpcdf=.false.                 !! true for netcdf runoff input

crunoffdir="./runoff/"
crunoffpre="roff____"           !! defaults runoff file name roff____yyyymmdd.one
crunoffsuf=".one"

crunoffcdf="NONE"
crofcdfvar="runoff"
syearin=0                       !! netcdf input file start date (set to 0 when not used)
smonin=0
sdayin=0
linterpcdf=.false.              !! true for netcdf input matrix

! * change
rewind(nsetfile)
read(nsetfile,ninput)
if(regionthis==1) then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 4. set_input****"
  write(nulnam,*) "linterp:"   ,linterp
  write(nulnam,*) "linterpcdf:"   ,linterpcdf
  write(nulnam,*) "cinpmat:",   trim(cinpmat)
  write(nulnam,*) "linpcdf:"   ,linpcdf
  write(nulnam,*) "crunoffdir:",trim(crunoffdir)
  write(nulnam,*) "crunoffpre:",trim(crunoffpre)
  write(nulnam,*) "crunoffsuf:",trim(crunoffsuf)
  write(nulnam,*) "crunoffcdf:",trim(crunoffcdf)
  write(nulnam,*) "crofcdfvar:",trim(crofcdfvar)
  write(nulnam,*) "syearin,smonin,sdayin",syearin,smonin,sdayin
endif

! ==============================
! *** 5. set output
! * defaults
coutdir="./"
crivoutdir="NONE"
crivstodir="NONE"
crivveldir="NONE"
crivdphdir="./"

cfldoutdir="NONE"
cfldstodir="NONE"
cflddphdir="NONE"
cfldfrcdir="./"
cfldaredir="./"

csfcelvdir="./"
coutflwdir="./"
cstorgedir="./"
coutinsdir="NONE"

cpthoutdir="NONE"
cpthflwdir="NONE"

loutvec=.false.
loutcdf=.false.
loutyyyy=.true.

! * change
rewind(nsetfile)
read(nsetfile,noutput)

if(regionthis==1) then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 5. set_output****"
  write(nulnam,*) "coutdir:",trim(coutdir)
  write(nulnam,*) "crivoutdir:",trim(crivoutdir)
  write(nulnam,*) "crivstodir:",trim(crivstodir)
  write(nulnam,*) "crivveldir:",trim(crivveldir)
  write(nulnam,*) "crivdphdir:",trim(crivdphdir)

  write(nulnam,*) "cfldoutdir:",trim(cfldoutdir)
  write(nulnam,*) "cfldstodir:",trim(cfldstodir)
  write(nulnam,*) "cflddphdir:",trim(cflddphdir)
  write(nulnam,*) "cfldfrcdir:",trim(cfldfrcdir)
  write(nulnam,*) "cfldaredir:",trim(cfldaredir)

  write(nulnam,*) "csfcelvdir:",trim(csfcelvdir)
  write(nulnam,*) "coutflwdir:",trim(coutflwdir)
  write(nulnam,*) "cstorgedir:",trim(cstorgedir)
  write(nulnam,*) "coutinsdir:",trim(coutinsdir)

  write(nulnam,*) "cpthoutdir:",trim(cpthoutdir)
  write(nulnam,*) "cpthflwdir:",trim(cpthflwdir)

  write(nulnam,*) "loutvec:",loutvec
  write(nulnam,*) "loutcdf:",loutcdf
  write(nulnam,*) "loutyyyy:",loutyyyy
endif

if( loutcdf ) loutvec=.false.

! ==============================
! *** 6. set conf: grid tstep configuration
! * defaults
nx=1440                                     !! 15 minute resolution
ny=720
nlfp=10                                     !! 10 floodplain layer
dt=24*60*60                                 !! dt = 1day (automatically set by adaptive time step)

nxin=360                                    !! 1 degree input
nyin=180
inpn=1                                      !! maximum number of input grids corresponding to one cama-flood grid
dtin=24*60*60                               !! default: /day -> /sec
drofunit=1.d-3                              !! default: mm   -> m3/m2

west=-180.d0                                !! west, east, north, south edges of the domain
east= 180.d0
north= 90.d0
south=-90.d0

lfld=.true.     ! floodplain scheme active     (.false. for no floodplain )
lkine=.false.   ! use local inertial equation  (.true.  for kinematic wave)
ladpstp=.true.  ! use adaptive time step       (.false. for no aptive time step)
lfldout=.true.  ! floodplan flow active        (.false. for no floodplain flow)
lpthout=.false. ! floodplan path flow          (.false. for no floodplain path flow)

ldam=.false.    ! dam operation                (.true.  for dam operation: not supported)
lmapend=.false. ! map data     in different endian (.false. for no endian conversion)
linpend=.false. ! input runoff in different endian (.false. for no endian conversion)

lleapyr=.true.  ! false to always use 365days/year

! * change
rewind(nsetfile)
read(nsetfile,nconf)

if( cdiminfo/="NONE" )then
  write(nulnam,*) " reading dimention info: ", trim(cdiminfo)
  open(tmpnam,file=cdiminfo,form='formatted')
  read(tmpnam,*) nx
  read(tmpnam,*) ny
  read(tmpnam,*) nlfp
  read(tmpnam,*) nxin
  read(tmpnam,*) nyin
  read(tmpnam,*) inpn
  read(tmpnam,*) 
  read(tmpnam,*) west
  read(tmpnam,*) east
  read(tmpnam,*) north
  read(tmpnam,*) south
  close(tmpnam)
endif

if(regionthis==1)then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 6. set_conf****"
  write(nulnam,*) "nx, ny,nlfp",nx,ny,nlfp
  write(nulnam,*) "dt",dt
  write(nulnam,*) "west,east,north,south",west,east,north,south
  write(nulnam,*) "nxin,nyin,dtin,drofunit",nxin,nyin,dtin,drofunit
  write(nulnam,*) "lfld",lfld
  write(nulnam,*) "lkine",lkine
  write(nulnam,*) "ladpstp",ladpstp
  write(nulnam,*) "lfldout",lfldout
  write(nulnam,*) "lpthout",lpthout
  write(nulnam,*) "lleapyr",lleapyr
endif

dt_def=dt                                   !! the default dt for adaptope time step mode

! *check
if ( .not.lfld .and. .not.lkine ) then
  write(nulnam,*) "lfld=.false. & lkine=.false."
  write(nulnam,*) "no floodplain run only available with kinematic wave (lkine=.true.)"
  write(nulnam,*) "stop"
  stop 9
endif

if ( lkine .and. ladpstp ) then
  write(nulnam,*) "lkine=.true. & ladpstp=.true."
  write(nulnam,*) "adaptive time step only available with local inertial equation (lkine=.false.)"
  write(nulnam,*) "stop"
  stop 9
endif

if ( lkine .and. lpthout ) then
  write(nulnam,*) "lkine=.true. & lpathout=.true."
  write(nulnam,*) "bifurcation channel flow only available with local inertial equation (lkine=.false.)"
  write(nulnam,*) "stop"
  stop 9
endif

if( lrestcdf .and. lpthout )then
  lstoonly=.true.
  write(nulnam,*) "lrestcdf=.true. & lpathout=.true."
  write(nulnam,*) "lstoonly set to .true."
endif

if( .not. lpthout )then
  cpthoutdir='NONE'
  cpthflwdir='NONE'
endif

if( loutcdf )then
  cpthflwdir='NONE'
endif

coutvardir(1)=trim(crivoutdir)
coutvardir(2)=trim(crivstodir)
coutvardir(3)=trim(crivdphdir)
coutvardir(4)=trim(crivveldir)

coutvardir(5)=trim(cfldoutdir)
coutvardir(6)=trim(cfldstodir)
coutvardir(7)=trim(cflddphdir)
coutvardir(8)=trim(cfldfrcdir)
coutvardir(9)=trim(cfldaredir)

coutvardir(10)=trim(csfcelvdir)
coutvardir(11)=trim(coutflwdir)
coutvardir(12)=trim(cstorgedir)
coutvardir(13)=trim(coutinsdir)

coutvardir(14)=trim(cpthoutdir)
coutvardir(15)=trim(cpthflwdir)


! ==============================
! *** 7. set param: parameters
! * defaults
pmanriv=0.03d0                              !! manning coefficient river
pmanfld=0.10d0                              !! manning coefficient floodplain
pgrv   =9.8d0                               !! gravity accerelation
pdstmth=10000.d0                            !! downstream distance at river mouth [m]
pcadp  =0.7d0                                !! cfl coefficient
pminslp=1.d-5                               !! minimum slope (kinematic wave)

! * change
rewind(nsetfile)
read(nsetfile,nparam)

if(regionthis==1)then
  write(nulnam,*) " "
  write(nulnam,*) "****init_inputnam: 7. parameters****"
  write(nulnam,*) "pmanriv, pmanfld", pmanriv, pmanfld
  write(nulnam,*) "pgrv   ",          pgrv
  write(nulnam,*) "pdstmth",          pdstmth
  write(nulnam,*) "pcadp  ",          pcadp
  write(nulnam,*) "pminslp",          pminslp
endif


! ==============================
! *** close file 
close(nsetfile)




end subroutine init_inputnam

end module init_inputnam_mod
