;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
path    = "/home/ys17-19/renql/project/2021variability"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
var_name = (/"preci"/)
unit     = (/"mm/day"/)
draw_var = (/"preci"/)      
ca = 80.0 ;day
cb = 8.0
month_s = 6
month_e = 8
nv = 0
siglvl = 0.05
sig_option  = False;True ;if False, only the significant values are shaded 
barshare = True ;False;
barvert  = True ;False;

if(case_option.eq.0) then
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
end if
pre_case   = (/"(a","(b","(c","(d","(a","(b","(g","(h","(i","(j","(k","(l"/)
fig_out    = "pdf";"x11";"eps";
fig_name = path + "/fig/" + case1(0) + "-corr_" + draw_var(nv) + "_" + cb + "-" + ca

;read preci data and then filter
lats1 = (/ 27, 32 , 25 , 40 /)
latn1 = (/ 37, 32 , 25 , 40 /)
lonl1 = (/105, 118, 115, 115/)
lonr1 = (/125, 118, 115, 115/)
na1 = 1

;area: east Asia
lats = (/0  ,0  ,10 /);-20
latn = (/90 ,70 ,55 /);60
lonl = (/0  ,0  ,90 /);20
lonr = (/360,180,140/);220
lat_spacing = 15
lon_spacing = 20
na = 2

f    = addfile(filein(3)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ntime = dimsizes(select_time)
ncase = dimsizes(case1)
nyear = dimsizes(years)
narea = dimsizes(lats1)

;the number of the panel plot
    plot   = new((/28/),graphic)
    plot2  = new((/28/),graphic)
    plot3  = new((/28/),graphic)
    wind   = new((/28/),graphic)
    wind2  = new((/28/),graphic)
    gtop   = new((/28/),graphic)
    line   = new((/28/),graphic)
    line1  = new((/28/),graphic)
    axis   = new((/28/),graphic)
    jetl   = new((/28/),graphic)
    sig    = new((/28/),graphic)

;set the area of the box
maxlon = (/lonr(na),lonr(na)/)
minlon = (/lonl(na),lonl(na)/)
mdlon1 = (/120,120/)
mdlon2 = (/240,240/)
maxlat = (/ 15,50/) 
minlat = (/-15,45/)
if(na.eq.1) then 
lonx = (/minlon,mdlon1,mdlon2,maxlon,maxlon,mdlon2,mdlon1,minlon,minlon/)
latx = (/minlat,minlat,minlat,minlat,maxlat,maxlat,maxlat,maxlat,minlat/)
else
lonx = (/minlon,maxlon,maxlon,minlon,minlon/)
latx = (/minlat,minlat,maxlat,maxlat,minlat/)
end if
printVarSummary(lonx)

;define the thickness and height of the title,tick and label
if(na.eq.0) then 
   tick_font      = 0.02
   subtitle_font  = 0.02
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

;=====================================================
;read data and filter and regression
;==================================================
wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file
do na1 = 0, narea-1, 1
np = 0
do nc = 0, ncase-1, 1
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

preci = new((/nyear,ntime/),float)
var := new((/nyear,ntime,nlat,nlon/),float)
if(nc.eq.0) then ;CN051
    f   = addfile(filein(nc),"r")
    t   = cd_calendar(f->time,0)
    select_time1 = ind(t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1) \
                    .and.t(:,1).ge.month_s.and.t(:,1).le.month_e)
    var = reshape(f->pre(select_time1,{lats(na):latn(na)},{lonl(na):lonr(na)}),(/nyear,ntime,nlat,nlon/));mm/day
    preci = reshape(wgt_areaave(f->pre(select_time1,{lats1(na1):latn1(na1)},{lonl1(na1):lonr1(na1)}),1.0,1.0,0),(/nyear,ntime/));mm/day
end if
if(nc.eq.1) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".nc","r")
    if(isleapyear(years(nt))) then ;leap year,366 day
    var(nt,:,:,:) = f->precip(select_time+1,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    preci(nt,:) = wgt_areaave(f->precip(select_time+1,{lats1(na1):latn1(na1)},{lonl1(na1):lonr1(na1)}),1.0,1.0,0)
    else
    var(nt,:,:,:) = f->precip(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    preci(nt,:) = wgt_areaave(f->precip(select_time,{lats1(na1):latn1(na1)},{lonl1(na1):lonr1(na1)}),1.0,1.0,0);mm/day
    end if
    end do
end if
if(nc.ge.2) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    var(nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    preci(nt,:) = wgt_areaave(f->PRECC(select_time,{lats1(na1):latn1(na1)},{lonl1(na1):lonr1(na1)}) ,1.0 ,1.0 ,0)
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    preci(nt,:) = preci(nt,:) + wgt_areaave(f->PRECL(select_time,{lats1(na1):latn1(na1)},{lonl1(na1):lonr1(na1)}) ,1.0 ,1.0 ,0)
    end do
    var = var*1000*3600*24 ;convert m/s to mm/day
    preci = preci*1000*3600*24 ;convert m/s to mm/day
end if
printVarSummary(var)
nsample = nyear*ntime
fca = 1.0/ca
fcb = 1.0/cb
opt = False
ts  = bw_bandpass_filter(preci,fca,fcb,opt,1)
var = bw_bandpass_filter(var,fca,fcb,opt,1)
corr := var(0,0,:,:)
corr = 0
do nt = 0, nyear-1, 1
    corr = corr + escorc_n(ts(nt,:), var(nt,:,:,:) ,0 ,0)
end do
corr = corr/nyear
printVarSummary(corr)
;ts := reshape(ts,(/nyear*ntime/))
;var:= reshape(var,(/nyear*ntime,nlat,nlon/))
;corr := escorc_n(ts, var ,0 ,0)
;prob := rtest(corr,nsample,0)
;if(sig_option) then 
;prob = mask(prob,prob.lt.siglvl,True)
;else
;corr = mask(corr,prob.lt.siglvl,True)
;end if
;copy_VarMeta(vars,prob)
copy_VarMeta(vars,corr)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
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

;----------Set res sig-----------
    res_sig                  = True       ; use plot options
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@gsnFrame         = False
    res_sig@gsnDraw          = False
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnLineThicknessF = 2
    ;res_sig@cnLineColor     = "black"
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.25
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

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
    resp@txString      = cb+"-"+ca+"-day filter, month: "+month_s+"-"+month_e+", preci (" + lats1(na1) + "-" + latn1(na1) + "N," + lonl1(na1) + "-" + lonr1(na1) + "E)"

minlvl = (/-0.8/)  ;q1 ;-80,-80
maxlvl = (/ 0.8/)      ; 80, 80
spacig = (/ 0.1/)      ; 10, 10

res@cnFillPalette = "BlueDarkRed18"
res@cnFillColors  = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)

    res@gsnLeftString   = case2(nc) 
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ") " ;
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np)    = gsn_csm_contour_map_ce(wks,corr,res)
    if(sig_option) then 
    sig(np) = gsn_csm_contour(wks,prob,res_sig)
    overlay(plot(np),sig(np))
    end if

    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    if(nc.ge.2) then 
    if(nc.eq.2.or.nc.eq.5) then
    line(np)  = gsn_add_polyline( wks, plot(np), lonx(:,0), latx(:,0), lnres )
    line1(np) = gsn_add_polyline( wks, plot(np), lonx(:,1), latx(:,1), lnres )
    else
    line(np) = gsn_add_polyline( wks, plot(np), lonx(:,nc-3), latx(:,nc-3), lnres )
    end if
    end if
    np = np + 1
    
    printMinMax(corr,True)
end do 
gsn_panel(wks,plot,(/3,2/),resp)
end do
end

