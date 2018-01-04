;******************************************************************************
;filename:      calc_clim_monthly_var.ncl
;output_file:   clim.mon.varname.nc
;calculate the global climate monthly variables
;                                            by Ql_Ren
;                                           2018/01/04
;******************************************************************************
begin
;------the path of the file
var_name = (/"omega"/)
filein   = "/home/ys17-19/renql/project/2018tej/data/NECP_DOE_monthly/" + var_name + ".mon.mean.nc"
fileout  = "/home/ys17-19/renql/project/2018tej/BaseMap/mdata/clim.mon." + var_name + ".nc"
years    = ispan(1979,2016,1)   ;the calculate time
nyear = dimsizes(years)

;---------------------------------------------------------------
;read the data, then calculate the climate monthly wind by clmMonTLLL
;--------------------------------------------------------------
f      = addfile(filein,"r")
var    = short2flt( f->$var_name(0)$(0:455,:,:,:) ) 
var_clim = clmMonTLLL(var)
printVarSummary(var_clim)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "Climate Monthly NCEP/DOE Reanalysis 2 from 1979 to 2016"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","level","lat","lon"/)
dimSizes = (/12,dimsizes(var&level),dimsizes(var&lat),dimsizes(var&lon)/)
dimUnlim = (/False, False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month","integer",(/"month"/))
filevardef(ncdf,"level",typeof(var&level),getvardims(var&level))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,var_name,typeof(var_clim),(/"month","level","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"level",var&level)
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,var_name,var_clim)

ncdf->month   = ispan(1,12,1) 
ncdf->level   = (/var&level/) 
ncdf->lat     = (/var&lat/)
ncdf->lon     = (/var&lon/)
ncdf->$var_name$    = (/var_clim/)
end

