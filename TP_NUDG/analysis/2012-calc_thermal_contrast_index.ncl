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
path    = "/home/ys17-19/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years    = ispan(1979,2005,1)   ;the calculate time
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
end if
if(case_option.eq.2) then
    case    = (/"F2000","NGTR"/) ;
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00"/)
                      ;"/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
end if
fileout = path + "/project/TP_NUDG/analysis/mdata/" + case(0) + "-thermal_contrast_index.txt"
var_name = (/"T","Q1"/)
draw_var = (/"uppT","intQ1"/)
nv = 1

;area: east Asia
area = (/"NH","TP","EA","NP"/)
lats = (/25  ,25  ,25  , 25 /);20 ;-20
latn = (/40  ,40  ,40  , 40 /) ;55 ;60
lonl = (/70  ,70  ,110 ,140 /) ;30 ;20
lonr = (/170 ,105 ,140 ,170 /) ;150;220
;lev  = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels ,hPa U
lev  = (/1000,925,850,700,600,  500,400,350,300,250, 200/) ;hPa
;lev  = (/1000,925,850,700,600,  500,400,350,300,250, 200,150,100/) ;hPa
;lev  = (/850,500,200/) ;v,

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev(0)},{lats(0):latn(0)},{lonl(0):lonr(0)})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
narea = dimsizes(area)
ncase = dimsizes(filein)
ntime = 365
var   = new((/ncase,nyear,365,nlat,nlon/),float)
g  = 9.8 ;m/(s*s)

do nf = 0,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
files2 := systemfunc("ls "+filein(nf)+"*.daily.PS.nc" )
do nt = 0,nyear-1,1
    print("handle with " + years(nt))
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f  = addfile(files(nt),"r")
    q1 = tofloat(f->$var_name(nv)$(:,{lev},{lats(0):latn(0)},{lonl(0):lonr(0)}))
    f  = addfile(files2(nt),"r")
    ps = f->PS(:,{lats(0):latn(0)},{lonl(0):lonr(0)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
    dp = dpres_plevel_Wrap(lev*100,ps,10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
    if(nv.eq.0) then 
    var(nf,nt,:,:,:) = dim_sum_n(q1(:,5:10,:,:)*dp(:,5:10,:,:),1)/dim_sum_n(dp(:,5:10,:,:),1)
    else
    var(nf,nt,:,:,:) = dim_sum_n(q1*dp,1)/g
    end if
end do ;nyear
end do ;ncase
delete([/q1,ps,dp,g/])
copy_VarMeta(vars,var(0,0,0,:,:))
var@long_name = "500-200hPa T"
var@units     = "K"
var_name(nv) = "uppT"
printVarSummary(var)

time = cd_calendar(f->time,0)
print("prob less than siglvl means significance")
do na = 2,narea-1,1
print("=================================================")
print(draw_var(nv)+" thermal contrast = "+area(1 )+"("+lats(1 )+"-"+latn(1 )+"N,"+lonl(1 )+"-"+lonr(1 )+"E) - " \
                                         +area(na)+"("+lats(na)+"-"+latn(na)+"N,"+lonl(na)+"-"+lonr(na)+"E)")
do time_scheme=0,2,1
if(time_scheme.eq.0) then 
months := (/"DJF","MAM","JJA","SON"/)
end if
if(time_scheme.eq.1) then 
months := (/"JF","MA","MJ","JA"/)
end if
if(time_scheme.eq.2) then 
months := (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
end if
nmonth = dimsizes(months) 
mons = 1
do nm = 0, nmonth-1, 1
    if(time_scheme.eq.0) then 
    if(nm.eq.0) then 
    month_s = 12  
    month_e = 2
    select_time := ind(time(:,1).eq.12.or.time(:,1).le.2)
    else
    month_s = 3*nm ;+ 1
    month_e = 3*nm + 2
    select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
    end if
    end if
    if(time_scheme.eq.1) then 
    month_s = 2*nm + mons
    month_e = 2*nm + mons+1
    select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
    end if
    if(time_scheme.eq.2) then 
    month_s = nm + mons
    month_e = nm + mons
    select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
    end if
    term  = dim_avg_n(var(:,:,select_time,:,:),2) ;get (nyear,nmonth,nlat,nlon)
    copy_VarMeta(vars,term(0,0,:,:))
    thcon = wgt_areaave(term(:,:,{lats(1):latn(1)},{lonl(1):lonr(1)}),1.0,1.0,0) - \
            wgt_areaave(term(:,:,{lats(na):latn(na)},{lonl(na):lonr(na)}),1.0,1.0,0) ;(ncase,nyear)
    delete(term)
    vari  = dim_variance_n(thcon,1)
    ave   = dim_avg_n(thcon,1)
    prob  = ttest(ave(0),vari(0),nyear,ave(1),vari(1),nyear,True,False)
    diff  = ave(1)-ave(0)
    print(months(nm)+": "+case(0)+": "+ave(0)+","+vari(0)+"  ; "+case(1)+": "+ave(1)+","+vari(1)+"  ; N-F diff: "+diff+"  ; prob: "+prob)
end do
print("")
end do
end do

end

