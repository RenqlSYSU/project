;***************************************************************************
;filename:      calc_clim_month_ave_wind.ncl
;output_file:   Clim_month_ave_wind.nc 
;1.read the 27 years variable(TLLL) of one case 
;2.calculate the month average to get an array of (nmonth,nlat,nlon),
;3.calculate the season average (nseason,nlat,nlon).
;
;Modification:
;20190712: add the case of noTP
;                                            by Ql_Ren
;                                           2018/01/18
;******************************************************************************
begin
path     = "/home/ys17-19"
;parameter for bw_bandpass_filter
ca = 80.0;90.0 ;day
cb = 8.0 ;2.5  ;10.0 ;
var_name = (/"U"  ,"V"  ,"T"  ,"Z3" ,"preci", \
             "preci-"+cb+"-"+ca+"day","TS" ,"PS" ,"PSL" /) ;,"OMEGA","RELHUM" 
case_option = 0
nv = 5
;area: east Asia
lats = 14 ;-90  ;0  
latn = 55 ; 90  ;60 
lonl = 69 ; 0   ;20 
lonr = 140; 360 ;180
lev  = (/1000,925,850,700,500,400, 300,200,100,50,10/) ;,150,100hPa
;lev  = (/850,500,200/) ;v,

if(case_option.eq.0) then
    case    = (/"CN05.1"/) ;
    filein  = path + "/data-observation/CN05.1/CN05.1_Pre_1961_2018_daily_025x025.nc"
    years   = ispan(1979,2005,1)   ;the calculate time
    f    = addfile(filein,"r")
    vars:= f->pre(0,{lats:latn},{lonl:lonr});mm/day
end if
if(case_option.eq.1) then
    case    = (/"APHRO"/) ;
    filein  = path + "/data-observation/APHRO_MA_025deg_V1003R1/APHRO_MA_025deg_V1003R1."; + years(nt) + ".nc" "
    years   = ispan(1979,2005,1)   ;the calculate time 1957~2007
    f    = addfile(filein+"2004.nc","r")
    vars:= f->precip(0,{lats:latn},{lonl:lonr});mm/day
    t    = cd_calendar(f->time,0)
    select_time2  = ind(.not.(t(:,1).eq.2.and.t(:,2).eq.29))
end if
fileout = path + "/renql/project/2021variability/data/" + case + "-Clim_month_ave_"

f1 = addfile("/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc","r")
dims  = dimsizes(vars)
nlat  = dims(0)
nlon  = dims(1)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
ncase = dimsizes(filein)
ntime = 365

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
if(nv.le.3) then 
    var        = new((/nyear,365,nlev,nlat,nlon/),float)
    var_month  = new((/nyear,12 ,nlev,nlat,nlon/),float)
    var_season = new((/nyear,4  ,nlev,nlat,nlon/),float)
    month_ave  = new((/      12 ,nlev,nlat,nlon/),float)
    month_vari = new((/      12 ,nlev,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
    season_ave = new((/       4 ,nlev,nlat,nlon/),float)
    season_vari= new((/       4 ,nlev,nlat,nlon/),float)
    
    time = cd_calendar(f1->time,0)
    do nm = 0,11,1
        select_time := ind(time(:,1).eq.nm+1)
        var_month(:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    end do
    ;----------calculate the season average-----------------------
    mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    do ns = 0,3,1
        nm = ns*3 ;nm = 0,3,6,9
        if(nm.eq.0) then
        var_season(:,ns,:,:,:) = ( var_month(:,11  ,:,:,:)*mdays(11)   + \
                                   var_month(:,nm  ,:,:,:)*mdays(nm)   + \
                                   var_month(:,nm+1,:,:,:)*mdays(nm+1))/(mdays(11) + mdays(nm)+ mdays(nm+1))
        else
        var_season(:,ns,:,:,:) = ( var_month(:,nm-1,:,:,:)*mdays(nm-1) + \
                                   var_month(:,nm  ,:,:,:)*mdays(nm)   + \
                                   var_month(:,nm+1,:,:,:)*mdays(nm+1))/(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
        end if
    end do
end if

if(nv.ge.4) then
    var        := new((/nyear,365,nlat,nlon/),float)
    var_month  := new((/nyear,12 ,nlat,nlon/),float)
    var_season := new((/nyear,4  ,nlat,nlon/),float)
    month_ave  := new((/      12 ,nlat,nlon/),float)
    month_vari := new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
    season_ave := new((/       4 ,nlat,nlon/),float)
    season_vari:= new((/       4 ,nlat,nlon/),float)
    
    if(case_option.eq.0) then
        f   = addfile(filein,"r")
        t   = cd_calendar(f->time,0)
        select_time1 = ind(.not.(t(:,1).eq.2.and.t(:,2).eq.29).and.t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1))
        print(dimsizes(select_time1))
        var = reshape(f->pre(select_time1,{lats:latn},{lonl:lonr}),(/nyear,365,nlat,nlon/));mm/day
    end if
    if(case_option.eq.1) then
        do nt = 0, nyear-1 ,1
        f   = addfile(filein+years(nt)+".nc","r")
        if(isleapyear(years(nt))) then ;leap year,366 day
        var(nt,:,:,:) = f->precip(select_time2,{lats:latn},{lonl:lonr});mm/day
        else
        var(nt,:,:,:) = f->precip(:,{lats:latn},{lonl:lonr});mm/day
        end if
        end do
    end if
    printVarSummary(var)
    print("missing value:" + num(ismissing(var)))
    var@_FillValue = -9999
    print("missing value:" + num(ismissing(var)))
    
    if(nv.eq.5) then 
    fca = 1.0/ca
    fcb = 1.0/cb
    opt = False
    var = var - conform(var,dim_avg_n(var,0),(/1,2,3/))  ;remove annual cycle
    var1 = linmsg_n(var,(/-1/),3)
    var1 = bw_bandpass_filter(var1,fca,fcb,opt,1)
    var1 = where(ismissing(var),var1@_FillValue,var1)
    var = var1*var1 
    var@long_name = "variance"
    var@units     = "mm2/day2" 
    end if

    time = cd_calendar(f1->time,0)
    do nm = 0,11,1
        select_time := ind(time(:,1).eq.nm+1)
        var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    end do
    ;----------calculate the season average-----------------------
    mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    do ns = 0,3,1
        nm = ns*3 ;nm = 0,3,6,9
        if(nm.eq.0) then
        var_season(:,ns,:,:) = ( var_month(:,11  ,:,:)*mdays(11)   + \
                                 var_month(:,nm  ,:,:)*mdays(nm)   + \
                                 var_month(:,nm+1,:,:)*mdays(nm+1))/(mdays(11) + mdays(nm)+ mdays(nm+1))
        else
        var_season(:,ns,:,:) = ( var_month(:,nm-1,:,:)*mdays(nm-1) + \
                                 var_month(:,nm  ,:,:)*mdays(nm)   + \
                                 var_month(:,nm+1,:,:)*mdays(nm+1))/(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
        end if
    end do
end if

month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
season_ave  = dim_avg_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
season_vari = dim_variance_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
printVarSummary(month_ave)
printVarSummary(season_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout+var_name(nv)+".nc")
ncdf = addfile(fileout+var_name(nv)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the monthly and seasonal average"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","season","lev","lat","lon"/)
dimSizes = (/  12   ,  4     ,nlev ,nlat ,nlon/)
dimUnlim = (/False  ,False   ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
if(case_option.eq.0) then 
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
else
filevardef(ncdf,"lat",typeof(vars&latitude),(/"lat"/))
filevardef(ncdf,"lon",typeof(vars&longitude),(/"lon"/))
end if
if((nv.le.3)) then 
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_ave",typeof(season_ave),(/"season","lev","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lev","lat","lon"/))
filevarattdef(ncdf,"lev",vars&level)
else
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave",typeof(season_ave),(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))
end if

;Define the attribute of the variables
if(case_option.eq.0) then 
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
else
filevarattdef(ncdf,"lat",vars&latitude)
filevarattdef(ncdf,"lon",vars&longitude)
end if
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"season_ave",var)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

ncdf->month  = ispan(1,12,1)
ncdf->season  = (/1212,345,678,91011/)
ncdf->lev = (/lev/)
if(case_option.eq.0) then 
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
else
ncdf->lat = (/vars&latitude/)
ncdf->lon = (/vars&longitude/)
end if
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end

