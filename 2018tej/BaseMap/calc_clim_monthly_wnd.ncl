;******************************************************************************
;filename:      calc_clim_monthly_wnd.ncl
;output_file:   clim.mon.uvwnd.nc
;calculate the global climate monthly wind of u and v
;                                            by Ql_Ren
;                                           2018/01/02
;******************************************************************************
begin
;------the path of the file
filein   = "/home/ys17-19/renql/project/2018tej/data/NECP_DOE_monthly/"
fileout  = "/home/ys17-19/renql/project/2018tej/BaseMap/mdata/clim.mon.uvwnd.nc"
var_name = (/"uwnd","vwnd"/)
years    = ispan(1979,2016,1)   ;the calculate time
nyear = dimsizes(years)

;---------------------------------------------------------------
;read the data, then calculate the climate monthly wind by clmMonTLLL
;--------------------------------------------------------------
f      = addfile(filein+var_name(0)+".mon.mean.nc","r")
var1   = short2flt( f->$var_name(0)$ ) 
var1_clim = clmMonTLLL(var1)
printVarSummary(var1_clim)

f      = addfile(filein+var_name(1)+".mon.mean.nc","r")
var2   = short2flt( f->$var_name(1)$ )
var2_clim = clmMonTLLL(var2)
printVarSummary(var2_clim)

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
dimSizes = (/12,dimsizes(var1&level),dimsizes(var1&lat),dimsizes(var1&lon)/)
dimUnlim = (/False, False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month","integer",(/"month"/))
filevardef(ncdf,"level",typeof(var1&level),getvardims(var1&level))
filevardef(ncdf,"lat",typeof(var1&lat),getvardims(var1&lat))
filevardef(ncdf,"lon",typeof(var1&lon),getvardims(var1&lon))
filevardef(ncdf,"uwnd",typeof(var1_clim),(/"month","level","lat","lon"/))
filevardef(ncdf,"vwnd",typeof(var2_clim),(/"month","level","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"level",var1&level)
filevarattdef(ncdf,"lat",var1&lat)
filevarattdef(ncdf,"lon",var1&lon)
filevarattdef(ncdf,"uwnd",var1_clim)
filevarattdef(ncdf,"vwnd",var2_clim)

ncdf->month   = ispan(1,12,1) 
ncdf->level   = (/var1&level/) 
ncdf->lat     = (/var1&lat/)
ncdf->lon     = (/var1&lon/)
ncdf->uwnd    = (/var1_clim/)
ncdf->vwnd    = (/var2_clim/)
end

