;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql/project/2021variability"
lats = (/ 27 , 22  , 30  /)
latn = (/ 37 , 30  , 45  /)
lonl = (/105 ,107  ,115  /)
lonr = (/125 ,130  ,130  /)
var_name = (/"preci","vort"/)
nv = 0 ;
lev = (/850,500,200,700,600,  100,10/)
nl = 1
stand = (/"","_stand"/)
ns = 0
allcase = True ;False ;if True, four case be draw in one figure
ifrunave= True ;False;
nave = 8 ;five days run ave
month_s = 5 
month_e = 9

first_areaave = True ;False;
siglvl = 0.1
na = 0 ;"EA","TP","SEA","NEA"
p_draw = (/120,90,60,40,30,20,15,10,7,5/)

if(nv.eq.0) then
    case1   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    years   = ispan(1979,2005,1)   ;the calculate time
    path2   = "/home/ys17-19"
    filein  = path2 +(/"/data-observation/CN05.1/CN05.1_Pre_1961_2018_daily_025x025.nc",\
                       "/data-observation/APHRO_MA_025deg_V1003R1/APHRO_MA_025deg_V1003R1.",\; + years(nt) + ".nc"
                       "/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
end if
if(nv.eq.1) then
    case1   = (/"ERAI","NCEP1","NCEP2","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"ERAI","NCEP1","NCEP2","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    years   = ispan(1979,2005,1)   ;the calculate time
    path2   = "/home/ys17-19"
    filein  = path2 +(/"/data-observation/ERA-interim/pressure/",\ ;/uwnd/uwnd.interim.1993.nc 1979-2018
                       "/data-observation/NCEP_NCAR/daily/",\ ;/uwnd/uwnd.1979.nc 1948-2018
                       "/data-observation/NCEP_DOE/daily/",\  ;/uwnd/uwnd.1979.nc,   short  1979-2020
                       "/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
end if
fig_name = path + "/fig/" + case1(0) + "-" + month_s + month_e + stand(ns) + "_area_spec_line_" + var_name(nv)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

f    = addfile(filein(3)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
ntime = dimsizes(select_time)
ncase = dimsizes(case1)
nyear = dimsizes(years)
nfrq = toint(ntime/2.0) 
indx = 4 ;the result of specx_ci

var  = new((/ncase,nyear,indx ,nfrq/),float)
ave  = new((/ncase      ,indx ,nfrq/),float)
plot = new((/ncase/),graphic)
vari2 = new((/ncase/),float)
vari3 = new((/ncase/),float)

;---------------------------------------------------------------
;read the data, then calculate the yearly power spectra 
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0.1;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

do nc = 0, ncase-1 ,1
print("")
print("handle with " + case1(nc))

;========= read data ======================================
if(nv.eq.0) then ;preci
if(nc.eq.0) then ;CN051
    f   = addfile(filein(nc),"r")
    vars:= f->pre(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
end if
if(nc.eq.1) then ;APHRO
    f   = addfile(filein(nc)+years(0)+".nc","r")
    vars:= f->precip(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
end if
if(nc.ge.2) then
    f   = addfile(filein(nc)+years(0)+".daily.PRECC.nc","r")
    vars:= f->PRECC(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
dim  = dimsizes(vars)
nlat = dim(0)
nlon = dim(1)

term := new((/nyear,ntime,nlat,nlon/),float)
if(nc.eq.0) then
    f   = addfile(filein(nc),"r")
    t   = cd_calendar(f->time,0)
    select_time1 = ind(t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1) \
                    .and.t(:,1).ge.month_s.and.t(:,1).le.month_e)
    term = reshape(f->pre(select_time1,{lats(na):latn(na)},{lonl(na):lonr(na)}),(/nyear,ntime,nlat,nlon/));mm/day
end if
if(nc.eq.1) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".nc","r")
    if(isleapyear(years(nt))) then ;leap year,366 day
    term(nt,:,:,:) = f->precip(select_time+1,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    else
    term(nt,:,:,:) = f->precip(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    end if
    end do
end if
if(nc.ge.2) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    term(nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    term(nt,:,:,:) = term(nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    end do
    term = term*1000*3600*24 ;convert m/s to mm/day
end if
end if

;-------------------vort------------------------------
if(nv.eq.1) then ;vort
if(nc.le.2) then 
if(nc.eq.0) then 
    f   = addfile(filein(nc)+"uwnd/uwnd.interim."+years(0)+".nc","r")
else
    f   = addfile(filein(nc)+"uwnd/uwnd."+years(0)+".nc","r")
end if
vars := f->uwnd(0,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.ge.3) then
    f   = addfile(filein(nc)+years(0)+".daily.PRECC.nc","r")
    vars:= f->PRECC(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
end if

dim  := dimsizes(vars)
nlat := dim(0)
nlon := dim(1)
term := new((/nyear,ntime,nlat,nlon/),float)
do nt = 0, nyear-1 ,1
if(nc.eq.0) then 
    f  = addfile(filein(nc)+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    u := f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    f  = addfile(filein(nc)+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    v := f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.eq.1) then 
    f  = addfile(filein(nc)+"uwnd/uwnd."+years(nt)+".nc","r")
    u := f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    f  = addfile(filein(nc)+"vwnd/vwnd."+years(nt)+".nc","r")
    v := f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.eq.2) then 
    f  = addfile(filein(nc)+"uwnd/uwnd."+years(nt)+".nc","r")
    u := short2flt(f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
    f  = addfile(filein(nc)+"vwnd/vwnd."+years(nt)+".nc","r")
    v := short2flt(f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
end if
if(nc.ge.3)
    f  = addfile(filein(nc)+years(nt)+".daily.U.nc","r")
    u := f->U(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    f  = addfile(filein(nc)+years(nt)+".daily.V.nc","r")
    v := f->V(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
term(nt,:,:,:) = uv2vr_cfd(u,v,v&lat,v&lon,2) 
end do
term = term*10^5 ;10^-5 s^-1
end if

;print(term(0,10,4:9,10))
;========= data processing and calc spectrum ======================================
term  = linmsg_n(term, -1, 3)
term  = linmsg_n(term, -1, 2)
if (any(ismissing(term))) then
    print("term has missing values")
end if
term = term - conform(term,dim_avg_n(term,0),(/1,2,3/))  ;remove annual cycle
opt  = 1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
dimh = 1
if(ifrunave) then 
term = runave_n_Wrap(term, nave, opt, dimh)
end if
printVarSummary(term)

if(first_areaave) then 
    term := wgt_areaave_Wrap(term,1.0,1.0,0) ;area average
    if(ns.eq.1) then
    term = dim_standardize_n(term,1,1) ;dim_standardize_n(x, opt, dims)
    end if
   
    spcavg = new ((/nfrq/),typeof(term))
    spcavg = 0.0
    r1zsum = 0.0
    do nt = 0, nyear-1, 1
        ; calculate spectrum
        sdof := specx_anal(term(nt,:),d,sm,pct) ;returns the degrees of freedom as a scalar.
        spcavg = spcavg + sdof@spcx                ; sum spc of each segment
        r1     = sdof@xlag1                        ; extract segment lag-1
        r1zsum = r1zsum  + 0.5*log((1+r1)/(1-r1)) ; sum the Fischer Z
    end do 
    r1z  = r1zsum/nyear                 ; average r1z
    r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1); transform back, this is the mean r1
    spcavg  = spcavg/nyear              ; average spectrum
    df      = 2.0*nyear   ; deg of freedom
    df@spcx = spcavg      ; assign the mean spc
    df@frq  = sdof@frq
    df@xlag1= r1          ; assign mean lag-1
    ave(nc,:,:) = specx_ci(df, siglvl, 1.0-siglvl) 

else

    ave1  = new((/indx,nfrq,nlat,nlon/),typeof(term))
    if(ns.eq.1) then
    term = dim_standardize_n(term,1,1) ;dim_standardize_n(x, opt, dims)
    end if

    spcavg = new ((/nfrq/),typeof(term))
    do ny = 0, nlat-1, 1
    do nx = 0, nlon-1, 1
    spcavg = 0.0
    r1zsum = 0.0
    do nt = 0, nyear-1, 1
        ; calculate spectrum
        sdof := specx_anal(term(nt,:,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
        spcavg = spcavg + sdof@spcx                ; sum spc of each segment
        r1     = sdof@xlag1                        ; extract segment lag-1
        r1zsum = r1zsum  + 0.5*log((1+r1)/(1-r1)) ; sum the Fischer Z
    end do 
    r1z  = r1zsum/nyear                 ; average r1z
    r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1); transform back, this is the mean r1
    spcavg  = spcavg/nyear              ; average spectrum
    df      = 2.0*nyear   ; deg of freedom
    df@spcx = spcavg      ; assign the mean spc
    df@frq  = sdof@frq
    df@xlag1= r1          ; assign mean lag-1
    ave1(:,:,ny,nx) = specx_ci(df, siglvl, 1.0-siglvl) 
    end do 
    end do 
    ave(nc,:,:) = wgt_areaave_Wrap(ave1,1.0,1.0,0) ;area average
    term := wgt_areaave_Wrap(term,1.0,1.0,0) ;area average
end if
vari2(nc) = dim_avg_n(dim_variance_n(term,1),0)
vari3(nc) = dim_variance(ndtooned(term))
end do ;end case

varia = ave*(1.0/ntime)
varia(:,:,0) = varia(:,:,0)/2.0
varia(:,:,nfrq-1) = varia(:,:,nfrq-1)/2.0
frq = sdof@frq
ave = ave*conform(ave,frq,2)

;calc total variance
vari1 = dim_sum_n(varia(:,0,:),1)
print(vari1)
print(vari2)
print(vari3)

print(sdof)
printVarSummary(ave)
ave!2   = "frq"
ave&frq = log(frq)
ave@long_name = "power*frequency"
ave@units     = "1/day"
;ave@long_name = "frequency"
;ave@units     = "1/day"

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.018
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 2.0
   font_index     = 22

wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file

    resxy     = True
    ;resxy@vpWidthF  = 0.45
    ;resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "Period (day)"          ; xaxis
    resxy@tiYAxisString = "Power*frq"           ; yaxis
    ;resxy@tiYAxisString = "Variance"           ; yaxis
    ;resxy@trYLog  = True
    resxy@trYMinF = 0
    resxy@trYMaxF = max(ave(:,0,:))
    resxy@trXMinF = log(min(frq))
    resxy@trXMaxF = log(max(frq)) 
    resxy@xyComputeXMin = True
    ;resxy@trXLog  = True
    resxy@tmXBMode    = "Explicit"
    resxy@tmXBValues  = log(1.0/p_draw)
    resxy@tmXBLabels  = p_draw
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case2 + " " + vari1
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@lgLabelFont            = font_index
    resxy@pmLegendWidthF   = 0.1
    resxy@pmLegendHeightF  = 0.25
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = - 0.28   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.28    ; move ref vector along X-axis
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font 
    resxy@tmXBLabelFontThicknessF  = tick_thick
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

    ;resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    ;resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    ;resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    ;resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    ;resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    ;resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont        = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa"
    if(nv.eq.1) then 
    resxy@gsnRightString   = lev(nl) + "hPa " + var_name(nv) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
    else
    resxy@gsnRightString   = var_name(nv) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
    end if

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    resp@txString      = "sm="+sm+" pct="+pct+" runave="+ifrunave+" nave="+nave

;ip   = ind(frq.le.0.2)    ; all indices for "short" periods
ip   = ind(frq.ge.0.01)    ; all indices for "short" periods

np = 0
if(allcase) then
if(nv.eq.0)
    resxy@xyDashPatterns    = (/1,2,0,0,0,0/)
    resxy@xyLineThicknesses = (/3.0,3.0,3.0,5.0,5.0,5.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","black","black","red","blue","grey"/)
else
    resxy@xyDashPatterns    = (/3,1,2,0,0,0,0/)
    resxy@xyLineThicknesses = (/3.0,3.0,3.0,3.0,5.0,5.0,5.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","black","black","black","red","blue","grey"/)
end if
    resxy@pmLegendDisplayMode = "Always"
    resxy@gsnLeftString   = pre_case(np) + ") " ;+ case2(nc)
    plot(np) = gsn_csm_xy(wks, ave&frq, ave(:,0,:), resxy)     
    np = np + 1
    gsn_panel(wks,plot,(/1,1/),resp)
else
    resxy@xyDashPatterns    = (/0,0,1,3/)
    resxy@xyLineThicknesses = (/2.0,3.0,3.0,3.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","red","red","red"/)
    resxy@pmLegendDisplayMode = "NoCreate"
    do nc = 0, ncase-1, 1
    resxy@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + vari1(nc)
    plot(np) = gsn_csm_xy(wks, ave&frq, ave(nc,:,:), resxy)     
    np = np + 1
    end do
    gsn_panel(wks,plot,(/1,ncase/),resp)
end if
end
