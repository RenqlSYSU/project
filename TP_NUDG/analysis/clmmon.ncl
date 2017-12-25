;***********************************************************
;calculate monthly climatology data from a daily climatology data
;******************************************************************
begin
    vars = (/"PRECC","PRECL"/)
    inname = "../TP-NUDG/pro/clim/AMIP_C5PM_TP_NUDG.cam.h1.ESM.clim.daily."
    outname = "./AMIP_C5PM_TP_NUDG.cam.h1.ESM.clim.month."
    ;vars = (/"T","Q","OMEGA","Z3","U","V"/)
    ;inname = "../AMIP-CTRL/pro/clim/AMIP_C5PM.cam.h1.ESM.clim."
    ;outname = "./AMIP_C5PM.cam.h1.ESM.clim.month."

do i = 0,dimsizes(vars)-1,1
   filein  = inname + vars(i) + ".nc" 
   fileout = outname + vars(i) + ".nc"   ;save data
  print(filein)

   f    =  addfile(filein, "r")
   var  =  f->$vars(i)$
   nlat = dimsizes(var&lat)
   nlon = dimsizes(var&lon)
   ;nlev = dimsizes(var&level)

   mdays = new((/12/),integer)
  ; clm_mon = new((/12,nlev,nlat,nlon/),float)
   clm_mon = new((/12,nlat,nlon/),float)
   mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   clm_mon = 0
   day2 = 0
   do m = 0,11
     day2 = day2+mdays(m)-1
     day1 = day2-mdays(m)+1
     do j = day1,day2,1
        ;clm_mon(m,:,:,:) = clm_mon(m,:,:,:) + var(j,:,:,:)
        clm_mon(m,:,:) = clm_mon(m,:,:) + var(j,:,:)
     end do
     ;clm_mon(m,:,:,:) = clm_mon(m,:,:,:)/mdays(m) 
     clm_mon(m,:,:) = clm_mon(m,:,:)/mdays(m) 
    end do

;Define Variable Properties
    ;copy_VarMeta(var(0,:,:,:), clm_mon(0,:,:,:))
    copy_VarMeta(var(0,:,:), clm_mon(0,:,:))
    clm_mon!0 = "month"
    clm_mon&month = ispan(1,12,1)

   ncdf = addfile(fileout,"c")
   ncdf->$vars(i)$ = clm_mon
   
end do

end 
