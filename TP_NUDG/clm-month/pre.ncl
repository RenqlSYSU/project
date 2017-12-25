begin
    vars = (/"PRECC","PRECL"/)
    inname = (/"./AMIP_C5PM.cam.h1.ESM.clim.month."\
        ,"./AMIP_C5PM_TP_NUDG.cam.h1.ESM.clim.month."/)
    ;vars = (/"T","Q","OMEGA","Z3","U","V"/)
    fileout = "/users/yangsong3/Model_Team/F/clm-month/preci.nc"

;*****************************************************
;define a new variable
;****************************************************
   f    =  addfile(inname(0)+vars(0)+".nc", "r")
   var  =  f->$vars(0)$
   nlat = dimsizes(var&lat)
   nlon = dimsizes(var&lon)
   ;nlev = dimsizes(var&level)
   ntime = dimsizes(var&month)
   preci = new((/dimsizes(inname),ntime,nlat,nlon/),float)

;***************************************************************************
;read the precc and precl from two file and add them
;***********************************************************************
preci = 0
do i = 0, dimsizes(inname)-1, 1
  do j = 0, dimsizes(vars)-1, 1
    filein = inname(i) + vars(j) + ".nc"
    f := addfile(filein,"r")
    preci(i,:,:,:) = preci(i,:,:,:) +  f->$vars(j)$
  end do
end do

preci = preci*1000*60*60*24 ;transform the unit from m/s to mm/day
copy_VarMeta(var,preci(0,:,:,:))
preci!0 = "case"
preci&case = (/"CTRL","NUDG"/)
print(preci(1,2,80,60))

system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@title = "daily precipitation in the climate month (mm/day)"
fileAtt@description = "data from two model case: 1 AMIP-CTRL ; 2 TP-NUDG."
fileattdef(ncdf,fileAtt)

dimNames = (/"month","lat","lon"/)
dimSizes = (/ntime,nlat,nlon/)
dimUnlim = (/False, False, False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

filevardef(ncdf,"month","integer","month")
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"CTRLpreci",typeof(var),getvardims(var))
filevardef(ncdf,"NUDGpreci",typeof(var),getvardims(var))
filevardef(ncdf,"DIFFpreci",typeof(var),getvardims(var))

filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"CTRLpreci",var)
filevarattdef(ncdf,"NUDGpreci",var)
filevarattdef(ncdf,"DIFFpreci",var)

ncdf->month = (/var&month/)
ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->CTRLpreci = (/preci(0,:,:,:)/)
ncdf->NUDGpreci = (/preci(1,:,:,:)/)
ncdf->DIFFpreci = (/preci(0,:,:,:)-preci(1,:,:,:)/)
end   

