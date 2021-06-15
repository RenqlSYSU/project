;******************************************************************************
;calculate the eof of observation and model daily summer precipitation, 
;which has been band filter
;use the function of eofunc_n_Wrap, eofunc_ts_n_Wrap to 
;calculate the EOF pattern and time serises.
;                                            by Ql_Ren
;                                           2021/06/08
;******************************************************************************
begin
path    = "/home/ys17-19/renql/project/2021variability"
var_name = (/"preci"/)
unit     = (/"mm/day"/)
nv = 0
ifrunave= 2 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
nave = 5 ;five days run ave
ca = 80.0 ;day
cb = 12.0
month_s = 6
month_e = 8
savedata = False;True
barshare = True ;False;
barvert  = True ;False;

opt = True
opt@jopt = 0 ;use covariance matrix (0) or correlation matrix (1) to compute EOFs
neval = 4
siglvl = 0.1

area = (/"EA"/);,"GB"
lats = (/ 20 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 45 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/105 /);, 0   ;30 ;20
lonr = (/135 /);, 360 ;150;220
lat_spacing = 10
lon_spacing = 15
na = 0

case1   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
case2   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
gtopfile= "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
years   = ispan(1979,2005,1)   ;the calculate time
path2   = "/home/ys17-19"
filein  = path2 +(/"/data-observation/CN05.1/CN05.1_Pre_1961_2018_daily_025x025.nc",\
                   "/data-observation/APHRO_MA_025deg_V1003R1/APHRO_MA_025deg_V1003R1.",\; + years(nt) + ".nc"
                   "/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                   "/renql/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                   "/renql/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                   "/renql/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
fileout  = path + "/data/eof_ts_"+area(na)+"_"+month_s+month_e+"_"+var_name(nv)+"_"+cb+"-"+ca+".nc"
fig_name = path + "/fig/" + case1(0) + "-eof_"+month_s+month_e+"_"+var_name(nv)+"_"+cb+"-"+ca
pre_case   = (/"(a","(b","(c","(d","(a","(b","(g","(h","(i","(j","(k","(l"/)
fig_out    = "pdf";"x11";"eps";

f    = addfile(filein(3)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
f = addfile(gtopfile,"r")
ldmk   = f->LANDM_COSLAT({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ntime = dimsizes(select_time)
ncase = dimsizes(case1)
nyear = dimsizes(years)
nsample = nyear*ntime
eof_ts  = new((/ncase,neval,nsample/),float)
ts_mean = new((/ncase,neval/),float)
plot    = new((/18/),graphic)
gtop    = new((/18/),graphic)

;define the thickness and height of the title,tick and label
if(na.eq.0) then 
   tick_font      = 0.02
   subtitle_font  = 0.025
   title_font     = 0.02
   vcref_font     = 0.008
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 1.5
   font_index     = 22
else
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.02
   vcref_font     = 0.015
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.5 ;3.0;
   vc_thick       = mp_thick + 0.5 
   font_index     = 22
end if
   label_font     = 0.018
   label_thick    = 1.5

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "landscape" ;"portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
    res@gsnAddCyclic     = False
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
   
;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "blueviolet";"black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = mp_thick + 1.0 
    ;lnres@gsLineLabelString    = "TP"

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = mp_thick + 1.5 
    res_gtop@cnLineColor      = "black"
    res_gtop@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_gtop@cnLevels         = (/3000/)
    res_gtop@gsnCenterString      = ""
    res_gtop@gsnLeftString    = " "
    res_gtop@gsnRightString   = " "

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = True
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.7
    else
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.8 ;0.05
    resp@pmLabelBarHeightF = 0.05;0.8
    end if
    else
    res@lbLabelBarOn         = True   ; the label bar
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarZone               = 0
    if(barvert) then
    res@lbOrientation        = "Vertical";
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.17
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.04
    res@pmLabelBarOrthogonalPosF  = 0.2 ;y direction
    res@pmLabelBarParallelPosF    = 0.5
    end if
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    if(ifrunave.eq.0) then 
    resp@txString      = "month: "+month_s+"-"+month_e+" no filter"
    end if
    if(ifrunave.eq.1) then 
    resp@txString      = "month: "+month_s+"-"+month_e+" runave="+nave+"day"
    end if
    if(ifrunave.eq.2) then 
    resp@txString      = "month: "+month_s+"-"+month_e+" bw_filter="+cb+"-"+ca+"day"
    end if
    if(ifrunave.eq.3) then 
    resp@txString      = "month: "+month_s+"-"+month_e+" Fourier_filter="+cb+"-"+ca+"day"
    end if

minlvl = (/-8/)  ;q1 ;-80,-80
maxlvl = (/ 8/)      ; 80, 80
spacig = (/ 1/)      ; 10, 10

res@cnFillPalette = "BlueDarkRed18"
res@cnFillColors  = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)

p_draw = (/120,90,60,40,30,20,15,10,7,5/)
    
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
    ;resxy@trYMinF = 0
    ;resxy@trYMaxF = max(ave(:,0,:))
    resxy@trXMinF = log(1.0/90.0) ;log(min(frq))
    resxy@trXMaxF = log(1.0/8.0)  ;log(max(frq)) 
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
    ;resxy@xyExplicitLegendLabels = case2 + " " + vari1
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
    resxy@gsnRightString   = var_name(nv) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"

; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0.1;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

;=====================================================
;read data and filter and regression
;==================================================
wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file
np = 0
do nc = 0, ncase-1, 1
if(nc.eq.0) then ;CN051
    f   = addfile(filein(nc),"r")
    vars:= f->pre(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    lat := vars&lat
end if
if(nc.eq.1) then ;APHRO
    f   = addfile(filein(nc)+years(0)+".nc","r")
    vars:= f->precip(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    lat := vars&latitude
end if
if(nc.ge.2) then
    f   = addfile(filein(nc)+years(0)+".daily.PRECC.nc","r")
    vars:= f->PRECC(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
    lat := vars&lat
end if
dim  = dimsizes(vars)
nlat = dim(0)
nlon = dim(1)

var := new((/nyear,ntime,nlat,nlon/),float)
if(nc.eq.0) then ;CN051
    f   = addfile(filein(nc),"r")
    t   = cd_calendar(f->time,0)
    select_time1 = ind(t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1) \
                    .and.t(:,1).ge.month_s.and.t(:,1).le.month_e)
    var = reshape(f->pre(select_time1,{lats(na):latn(na)},{lonl(na):lonr(na)}),(/nyear,ntime,nlat,nlon/));mm/day
end if
if(nc.eq.1) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".nc","r")
    if(isleapyear(years(nt))) then ;leap year,366 day
    var(nt,:,:,:) = f->precip(select_time+1,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    else
    var(nt,:,:,:) = f->precip(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    end if
    end do
end if
if(nc.ge.2) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    var(nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    end do
    var = var*1000*3600*24 ;convert m/s to mm/day
end if
printVarSummary(var)

dimh = 1
if(ifrunave.eq.1) then 
opt1= -1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
var = runave_n(var, nave, opt1, dimh)
var = var - conform(var,dim_avg_n(var,0),(/1,2,3/))  ;remove annual cycle
end if
if(ifrunave.eq.2) then 
fca = 1.0/ca
fcb = 1.0/cb
opt2= False
var = bw_bandpass_filter(var,fca,fcb,opt2,dimh)
end if
pi   = atan(1.0)*4
wgt := sqrt(cos( lat*pi/180.0 ))
var = tofloat(var*conform(var,wgt,2))
var:= reshape(var,(/nyear*ntime,nlat,nlon/))

if(nc.ge.2) then 
    var = mask(var, conform(var,ldmk,(/1,2/)).eq.1, True)
end if
eof     := eofunc_n(var,neval,opt,0) ;last one is the dim ; get the (neval,nlat,nlon)
sig     := eofunc_north(eof@pcvar,nsample,True)
eof_ts_temp   = eofunc_ts_n(var,eof,False,0) ;get the (neval,nyear)
eof_ts(nc,:,:)= eof_ts_temp/conform(eof_ts_temp,sqrt(eof@eval),0)
ts_mean(nc,:) = eof_ts_temp@ts_mean
eof = eof*conform(eof,sqrt(eof@eval),0)
print(case1(nc)+" eval : " + eof@eval(0) + ", " + eof@eval(1) \
        + ", " + eof@eval(2) + ", " + eof@eval(3)) 
subtitle_color = where(sig,"red","black") 
copy_VarMeta(vars,eof(0,:,:))

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
do ne = 0, 2, 1
    nee = ne+1
    res@gsnLeftString   = case2(nc);+" " +var_name(nv) + " (" + unit(nv) + ") " ; 
    res@gsnRightStringFontColor   = subtitle_color(ne)
    res@gsnRightString            = "EOF"+nee+": "+eof@eval(ne)+"  "+eof@pcvar(ne)+"%" 
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np) = gsn_csm_contour_map_ce(wks,eof(ne,:,:),res)
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    np = np + 1
    printMinMax(eof(ne,:,:),True)
end do 
end do 
gsn_panel(wks,plot,(/6,3/),resp)

indx = 4 ;the result of specx_ci
nsample2 = ntime;nyear*
nfrq  = toint(nsample2/2.0) 
splt = new((/ncase,indx,nfrq/),float)
eofts = reshape(eof_ts,(/ncase,neval,nyear,ntime/))

np = 0
do ne = 0, 2, 1
nee = ne+1
do nc = 0, ncase-1, 1
if(nsample2.eq.ntime) then 
    spcavg = new ((/nfrq/),typeof(eofts))
    spcavg = 0.0
    r1zsum = 0.0
    do nt = 0, nyear-1, 1
        ; calculate spectrum
        sdof := specx_anal(eofts(nc,ne,nt,:),d,sm,pct) ;returns the degrees of freedom as a scalar.
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
    splt(nc,:,:) = specx_ci(df, siglvl, 1.0-siglvl) 
else
    sdof := specx_anal(eof_ts(nc,ne,:),d,sm,pct) ;returns the degrees of freedom as a scalar.
    splt(nc,:,:) = specx_ci(sdof, siglvl, 1.0-siglvl)
end if
end do
varia = splt*(1.0/nsample2)
varia(:,:,0) = varia(:,:,0)/2.0
varia(:,:,nfrq-1) = varia(:,:,nfrq-1)/2.0
frq  = sdof@frq
splt = splt*conform(splt,frq,2)
vari1 = dim_sum_n(varia(:,0,:),1)

print(sdof)
printVarSummary(splt)
splt!2   = "frq"
splt&frq = log(frq)
splt@long_name = "power*frequency"
splt@units     = "1/day"

    resxy@xyExplicitLegendLabels = case2 + " " + vari1
    resxy@xyDashPatterns    = (/1,2,0,0,0,0/)
    resxy@xyLineThicknesses = (/3.0,3.0,3.0,5.0,5.0,5.0/)-1.5   ; make second line thicker
    resxy@xyLineColors      = (/"black","black","black","red","blue","grey"/)
    resxy@pmLegendDisplayMode = "Always"
    resxy@gsnLeftString   = pre_case(np) + ") EOF"+nee
    plot(np) = gsn_csm_xy(wks, splt&frq, splt(:,0,:), resxy)     
    np = np + 1
end do
resp@gsnPanelLabelBar = False 
gsn_panel(wks,plot,(/1,3/),resp)

;========================================================================
;save the data
;=======================================================================
if(savedata) then 
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "eof of ISV preci"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"case","neval","time"/)
dimSizes = (/ ncase ,neval  ,nsample /)
dimUnlim = (/False  ,False  ,False /)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"case","character",(/"case"/))
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"time","integer",(/"time"/))
filevardef(ncdf,"eof_ts",typeof(eof_ts),(/"case","neval","time"/))

;Define var Attribute
filevarattdef(ncdf,"time",time)
filevarattdef(ncdf,"eof_ts",eof_ts)

ncdf->case  = case1 
ncdf->time  = ispan(1,nsample-1,1) 
ncdf->eof_ts = (/eof_ts/)
end if
end

