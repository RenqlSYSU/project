;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
path    = "/home/ys17-19/renql/project/2021variability"
var_name = (/"vort","OLR","preci"/)
unit     = (/"10~S~-6~N~ m~S~-1~N~"," W m~S~-2~N~","mm/day"/)
draw_var = (/"vort","OLR","preci"/)      
month_s = 5
month_e = 9
lev = (/850,500,200,700,600,  100,10/)
nl = 0
nv = 1
nlag = 6
siglvl = 0.05
sig_option  = False;True ;if False, only the significant values are shaded 
vect_option = True ;False;
barshare = True ;False;
barvert  = False;True ;

gtopfile= "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
years   = ispan(1979,2005,1)   ;the calculate time
path2   = "/home/ys17-19/data-observation"
case1   = (/"CN05.1","APHRO"/) 
filein1 = path2 +(/"/CN05.1/CN05.1_Pre_1961_2018_daily_025x025.nc",\
                   "/APHRO_MA_025deg_V1003R1/APHRO_MA_025deg_V1003R1."/); + years(nt) + ".nc" ;1951-2007
nc1 = 0
case2   = (/"ERAI","NCEP1","NCEP2"/)
filein2 = path2 +(/"/ERA-interim/pressure/",\ ;/uwnd/uwnd.interim.1993.nc 1979-2018
                   "/NCEP_NCAR/daily/",\ ;/uwnd/uwnd.1979.nc 1948-2018
                   "/NCEP_DOE/daily/"/)  ;/uwnd/uwnd.1979.nc,   short  1979-2020
pre_case   = (/"(a","(b","(c","(d","(a","(b","(g","(h","(i","(j","(k","(l"/)
fig_out    = "pdf";"x11";"eps";
fig_name = path + "/fig/" + case1(0) + "-obs_regr_" + lev(nl) + "wind_" + draw_var(nv)

;read preci data and then filter
lats1 = (/ 27 /);-20
latn1 = (/ 37 /);60
lonl1 = (/105 /);20
lonr1 = (/125 /);220

;area: east Asia
lats = (/0  ,0  /);-20
latn = (/90 ,70 /);60
lonl = (/0  ,0  /);20
lonr = (/360,180/);220
lat_spacing = 20
lon_spacing = 20
na = 1

f    = addfile("/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.PRECC.nc","r")
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
ncase1 = dimsizes(case1)
ncase2 = dimsizes(case2)
nyear = dimsizes(years)
nvar  = 3
preci = new((/nyear,ntime/),float)

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

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)

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
   
;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow";"CurlyVector";
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
if(na.eq.0) then 
    resvc@vcMinDistanceF            = 0.009
    resvc@vcRefLengthF              = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
else
    resvc@vcMinDistanceF            = 0.02
    resvc@vcRefLengthF              = 0.05 ;0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.008;0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.008;0.015
end if
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    resvc@vcRefAnnoOn                = False
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.3          ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False;True
    resvc@vcRefAnnoString2On         = False;True ;
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor           = "grey"

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "blueviolet";"black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = mp_thick + 1.0 
    ;lnres@gsLineLabelString    = "TP"

    res_axis                      = True
    res_axis@gsLineColor          = "blueviolet";"black"
    res_axis@gsLineDashPattern    = 1
    res_axis@gsLineThicknessF     = mp_thick + 2.5 
    ;res_axis@gsLineLabelString    = "Westly"

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

if(na.eq.0) then ; 850,700,600,500,200,100,10
    vcref         = (/5  ,5 ,5 ,5 /)  ;the level of 850hpa    ctrl,nudg,diff
else                               
    vcref         = (/0.5, 1, 1,1 /)  ;the level of 850hpa    ctrl,nudg,diff
end if

minlvl = (/-4.0,-1.6/)  ;q1 ;-80,-80
maxlvl = (/ 4.0, 1.6/)      ; 80, 80
spacig = (/ 0.5, 0.2/)      ; 10, 10

res@cnFillPalette = "BlueDarkRed18"
res@cnFillColors  = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)

;=====================================================
;read data and filter and regression
;==================================================
if(nc1.eq.0) then 
    f   = addfile(filein1(nc1),"r")
    t   = cd_calendar(f->time,0)
    select_time1 = ind(t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1) \
                    .and.t(:,1).ge.month_s.and.t(:,1).le.month_e)
    term  = f->pre(select_time1,{lats1:latn1},{lonl1:lonr1})
    dim  := dimsizes(term)
    nlat := dim(1)
    nlon := dim(2)
    term  := reshape(term,(/nyear,ntime,nlat,nlon/));mm/day
    preci = wgt_areaave_Wrap(term, 1.0, 1.0, 0)

    if(nv.eq.1) then 
    term := f->pre(select_time1,{lats(na):latn(na)},{lonl(na):lonr(na)})
    dim  := dimsizes(term)
    nlat1 := dim(1)
    nlon1 := dim(2)
    rain  = reshape(term,(/nyear,ntime,nlat1,nlon1/));mm/day
    end if
end if
if(nc1.eq.1) then 
    do nt = 0, nyear-1 ,1
    f   = addfile(filein1(nc1)+years(nt)+".nc","r")
    if(isleapyear(years(nt))) then ;leap year,366 day
    preci(nt,:) = wgt_areaave_Wrap(f->precip(select_time+1,{lats1:latn1},{lonl1:lonr1}),1.0,1.0,0);mm/day
    else
    preci(nt,:) = wgt_areaave_Wrap(f->precip(select_time,{lats1:latn1},{lonl1:lonr1}),1.0,1.0,0);mm/day
    end if
    end do
end if

do nc = 0, ncase2-1, 1
resp@txString      = case2(nc)+" month: "+month_s+"-"+month_e+", "+case1(nc1)+"preci ("+lats1+"-"+latn1+"N,"+lonl1+"-"+lonr1+"E)"
if(nc.eq.0) then 
    f   = addfile(filein2(nc)+"uwnd/uwnd.interim."+years(0)+".nc","r")
else
    f   = addfile(filein2(nc)+"uwnd/uwnd."+years(0)+".nc","r")
end if
vars := f->uwnd(0,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
dim  := dimsizes(vars)
nlat := dim(0)
nlon := dim(1)
var1 := new((/nvar,nyear,ntime,nlat,nlon/),float)

do nt = 0, nyear-1 ,1
if(isleapyear(years(nt))) then ;leap year,366 day
if(nc.eq.0) then 
    inv = 0
    f   = addfile(filein2(nc)+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->uwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc)+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->vwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.eq.1) then 
    inv = 0
    f   = addfile(filein2(nc)+"uwnd/uwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->uwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc)+"vwnd/vwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->vwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.eq.2) then 
    inv = 0
    f   = addfile(filein2(nc)+"uwnd/uwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = short2flt(f->uwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
    inv = 1
    f   = addfile(filein2(nc)+"vwnd/vwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = short2flt(f->vwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
end if
else
if(nc.eq.0) then 
    inv = 0
    f   = addfile(filein2(nc)+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc)+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.eq.1) then 
    inv = 0
    f   = addfile(filein2(nc)+"uwnd/uwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc)+"vwnd/vwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc.eq.2) then 
    inv = 0
    f   = addfile(filein2(nc)+"uwnd/uwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = short2flt(f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
    inv = 1
    f   = addfile(filein2(nc)+"vwnd/vwnd."+years(nt)+".nc","r")
    var1(inv,nt,:,:,:) = short2flt(f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
end if
end if
end do
if(nv.eq.0) then ;vort
    var1(2,:,:,:,:) = uv2vr_cfd(var1(0,:,:,:,:),var1(1,:,:,:,:),var1&lat,var1&lon,2) 
    var1(2,:,:,:,:) = var1(2,:,:,:,:)*10^6
end if

do lgop = 0, 3, 1
lagg = (/5,4,3,2,1,0/)
if(lgop.eq.0) then 
ca = 90.0 ;day
cb = 8.0
lagg = lagg
end if
if(lgop.eq.1) then 
ca = 20.0 ;day
cb = 15.0
lagg = 2*lagg
end if
if(lgop.eq.2) then 
ca = 25.0 ;day
cb = 17.0
lagg = 2*lagg
end if
if(lgop.eq.3) then 
ca = 60.0 ;day
cb = 25.0
lagg = 4*lagg
end if
fca = 1.0/ca
fcb = 1.0/cb
opt = False
;ts = bw_bandpass_filter(reshape(preci,(/ncase,nyear*ntime/)),fca,fcb,opt,1)
ts = reshape(bw_bandpass_filter(preci,fca,fcb,opt,1),(/nyear*ntime/))
var := reshape(bw_bandpass_filter(var1,fca,fcb,opt,2),(/nvar,nyear*ntime,nlat,nlon/))
;var := reshape(var1,(/ncase,nvar,nyear*ntime,nlat,nlon/))

if(nv.eq.1) then 
rain1 := reshape(bw_bandpass_filter(rain,fca,fcb,opt,2),(/nyear*ntime,nlat1,nlon1/))
regr1 := new((/nlag,nlat1,nlon1/),float)
tval1 := new((/nlag,nlat1,nlon1/),float)
df1   := new((/nlag,nlat1,nlon1/),float)
prob1 := new((/nlag,nlat1,nlon1/),float)
do nlg = 0, nlag-1, 1
    regr1(nlg,:,:) = regCoef_n(ts(lagg(nlg):(nyear*ntime-1)),rain1(0:(nyear*ntime-lagg(nlg)-1),:,:),0,0)
    tval1(nlg,:,:) = onedtond(regr1@tval ,(/nlat1,nlon1/))
    df1(nlg,:,:)   = onedtond(regr1@nptxy,(/nlat1,nlon1/))
end do
b := tval1 
b = 0.5
prob1  := betainc( df1/(df1+tval1^2), df1/2.0, b )
delete([/tval1,df1/])

regr1 = mask(regr1,prob1.lt.siglvl,True)
copy_VarCoords(term(0,:,:), regr1(0,:,:))
copy_VarCoords(term(0,:,:), prob1(0,:,:))
printVarSummary(regr1)
end if

regr := new((/nvar,nlag,nlat,nlon/),float)
tval := new((/nvar,nlag,nlat,nlon/),float)
df   := new((/nvar,nlag,nlat,nlon/),float)
prob := new((/nvar,nlag,nlat,nlon/),float)
do nlg = 0, nlag-1, 1
    regr(:,nlg,:,:) = regCoef_n(ts(lagg(nlg):(nyear*ntime-1)),var(:,0:(nyear*ntime-lagg(nlg)-1),:,:),0,1)
    tval(:,nlg,:,:) = onedtond(regr@tval ,(/nvar,nlat,nlon/))
    df(:,nlg,:,:)   = onedtond(regr@nptxy,(/nvar,nlat,nlon/))
end do
copy_VarCoords(vars, regr(0,0,:,:))
b := tval 
b = 0.5
prob  := betainc( df/(df+tval^2), df/2.0, b )
delete([/tval,df/])

notsiguv := mask(regr(0:1,:,:,:),prob(0:1,:,:,:).lt.siglvl,False)
if(sig_option) then
regr(0:1,:,:,:) = mask(regr(0:1,:,:,:),prob(0:1,:,:,:).lt.siglvl,True)
else
regr = mask(regr,prob.lt.siglvl,True)
end if
prob = mask(prob,prob.lt.siglvl,True)
copy_VarCoords(vars, notsiguv(0,0,:,:))
printVarSummary(regr)
printVarSummary(notsiguv)

do nlg = 0,nlag-1,1
    np = nlg*4+lgop
    if(lagg(nlg).eq.0) then 
    res@gsnLeftString   = cb+"-"+ca+"-day (" + lagg(nlg) + " day)" ;pre_case(np) + ") " +  
    else                                                                   
    res@gsnLeftString   = cb+"-"+ca+"-day (-" + lagg(nlg) + " day)";pre_case(np) + ") " +  
    end if
    res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ") " ;
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    if(nv.eq.1) then 
    plot(np)    = gsn_csm_contour_map_ce(wks,regr1(nlg,:,:),res)
    else
    plot(np)    = gsn_csm_contour_map_ce(wks,regr(2,nlg,:,:),res)
    end if
    if(sig_option) then 
    sig(np) = gsn_csm_contour(wks,prob(2,nlg,:,:),res_sig)
    overlay(plot(np),sig(np))
    end if

    if(vect_option) then
    resvc@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl)+"m/s"
    resvc2@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc2@vcRefAnnoString2     = vcref(nl)+"m/s"
    wind2(np)   = gsn_csm_vector(wks,notsiguv(0,nlg,:,:),notsiguv(1,nlg,:,:),resvc2)
    overlay(plot(np),wind2(np))
    wind(np)    = gsn_csm_vector(wks,regr(0,nlg,:,:),regr(1,nlg,:,:),resvc)
    overlay(plot(np),wind(np))
    end if
    
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    line(np)  = gsn_add_polyline( wks, plot(np), lonx(:,0), latx(:,0), lnres )
    line1(np) = gsn_add_polyline( wks, plot(np), lonx(:,1), latx(:,1), lnres )
    
    print(max(regr(2,nlg,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    print(min(regr(2,nlg,{lats(na):latn(na)},{lonl(na):lonr(na)})))
end do 
end do 
gsn_panel(wks,plot,(/6,4/),resp)
end do
end

