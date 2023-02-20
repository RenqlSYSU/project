;***************************************************************************
;convert F2000_CAM5.cam.h1.ESM.clim.T.nc from hybrid sigma pressure coordinate
;to pressure coordinate
;                                            by Ql_Ren
;                                           2020/11/28
;******************************************************************************
begin
;the path of the file
var_name= (/"U","V","T"/)
path    = "/WORK/sysu_hjkx_ys/renql/NUDG/"
filein  = path + "F2000_CAM5.cam.h1.ESM.clim." + var_name + ".nc"
fileout = path + "F2000-clim_" + var_name + ".nc" 
nvar = dimsizes(var_name)

fcoff = addfile(path+"F2000_CAM5_TP_NUDG/pre_data/F2000_CAM5_TP_NUDG.cam.h1.0011-11-12-00000.nc")
    hbcofa	    =   f->hyam(0:layers-1)
    hbcofb	    =   f->hybm(0:layers-1)
    ctrl_psfc   =   f->PS 
    p0	        =   1000

;---------------------------------------------------------------
;read the data, then convert 
;--------------------------------------------------------------
do nv = 0, nvar-1, 1 
    f    = addfile(filein(nv))
    data = f->$var_name(nv)$

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc))
ncdf = addfile(fileout(nc),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The annual cycle of TP average" 
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"year","day","lev"/)
dimSizes = (/ nyear,365  ,nlev /)
dimUnlim = (/ False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"day" ,"integer",(/"day"/))
filevardef(ncdf,"lev" ,"integer",(/"lev"/))
filevardef(ncdf,"series",typeof(var),(/"year","day","lev"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",var&level)
filevarattdef(ncdf,"series",var)
filevarattdef(ncdf,"ave",ave)

ncdf->day    = ispan(1,365,1)
ncdf->lev   = (/lev/)
ncdf->series = (/var/)
ncdf->ave = (/ave/)
end do
end


