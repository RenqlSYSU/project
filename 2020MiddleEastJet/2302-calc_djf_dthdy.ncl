;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path1   = "/home/ys17-23/Extension/ERA-interim/monthly/air.mon.mean.nc"
path2   = "/home/ys17-23/Extension2/renql/project/2020MiddleEastJet/data/"
var_name = (/"dthdy"/)

;area: east Asia
lats = 0  
latn = 90 
lonl = 0 
lonr = 360
lev  = (/925,850,500,200/)

years    = ispan(1979,2016,1)
years@calendar = "standard"
DJF_option = False

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f = addfile(path1,"r")
time = cd_calendar(f->time,0)
printVarSummary(time)

select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
printVarSummary(select_time)
vars  = f->t(select_time(2:112),{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&latitude)
nlon  = dimsizes(vars&longitude)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
var = dim_avg_n_Wrap(reshape(vars,(/nyear-1,3,nlev,nlat,nlon/)),1)

fileout  = path2 + "ERA-Interim_DJF_" + var_name + "-year.nc"

;------calc dthdy--------------------------------
    pi = atan(1.0)*4
    lat = (vars&latitude)*pi/180.0
    lev_dim = 1
    th  = pot_temp(lev*100, var, lev_dim, False)
    
    a  = 6378388 ;the radius of earth, m
    opt    = 0     ;used by center_finite_diff_n, no meanging 
    var = center_finite_diff_n(th,lat,False,opt,2)/a 
    printVarSummary(var)
    printMinMax(var,True)
    delete([/th,lats/])

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "DJF mean" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"year" ,"lev","lat","lon"/);
dimSizes = (/nyear-1,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&latitude),(/"lat"/))
filevardef(ncdf,"lon" ,typeof(vars&longitude),(/"lon"/))
filevardef(ncdf,"season_ave",typeof(var),(/"year","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&latitude)
filevarattdef(ncdf,"lon",vars&longitude)
filevarattdef(ncdf,"season_ave",var)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&latitude/)
ncdf->lon  = (/vars&longitude/)
ncdf->season_ave  = (/var/)
end

