SUBROUTINE init_diag

use parkind1 ,only: jpim, jprb, jprm
use mod_map  ,only: nseqmax, npthout, npthlev
use mod_diag ,only: d2rivinf, d2rivdph, d2rivvel, &
                    d2fldinf, d2flddph, d2fldfrc, d2fldare,&
                    d2pthout, d2pthinf, d2sfcelv, d2outflw, d2storge
use mod_diag ,only: d2outflw_avg, d2rivout_avg, d2fldout_avg, d2rivvel_avg, &
                    d2pthout_avg, d1pthflw_avg, d2outins, nadd

! ===============================================

allocate(d2rivinf(nseqmax,1))
allocate(d2rivdph(nseqmax,1))
allocate(d2rivvel(nseqmax,1))

allocate(d2fldinf(nseqmax,1))
allocate(d2flddph(nseqmax,1))
allocate(d2fldfrc(nseqmax,1))
allocate(d2fldare(nseqmax,1))

allocate(d2pthout(nseqmax,1))
allocate(d2pthinf(nseqmax,1))

allocate(d2sfcelv(nseqmax,1))
allocate(d2outflw(nseqmax,1))
allocate(d2storge(nseqmax,1))

allocate(d2rivout_avg(nseqmax,1))
allocate(d2fldout_avg(nseqmax,1))
allocate(d2rivvel_avg(nseqmax,1))
allocate(d2pthout_avg(nseqmax,1))
allocate(d2outflw_avg(nseqmax,1))    !! bugfix v362
allocate(d1pthflw_avg(npthout,npthlev))
allocate(d2outins(nseqmax,1))

d2rivinf(:,:)=0.d0
d2rivdph(:,:)=0.d0
d2rivvel(:,:)=0.d0

d2fldinf(:,:)=0.d0
d2flddph(:,:)=0.d0
d2fldfrc(:,:)=0.d0
d2fldare(:,:)=0.d0

d2pthout(:,:)=0.d0
d2pthinf(:,:)=0.d0

d2sfcelv(:,:)=0.d0
d2outflw(:,:)=0.d0
d2storge(:,:)=0.d0

d2rivout_avg(:,:)=0.d0
d2fldout_avg(:,:)=0.d0
d2rivvel_avg(:,:)=0.d0
d2pthout_avg(:,:)=0.d0
d1pthflw_avg(:,:)=0.d0
d2outflw_avg(:,:)=0.d0     !! bugfix v362
nadd=0

d2outins(:,:)=0.d0

END SUBROUTINE init_diag
