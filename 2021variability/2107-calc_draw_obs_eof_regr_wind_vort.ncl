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
unit     = (/"10~S~-6~N~ m~S~-1~N~"," W/m~S~2~N~","mm/day"/)
draw_var = (/"vort","OLR","preci"/)      
lev = (/850,500,200,700,600,  100,10/)
nl = 2
nv = 1

ifrunave2= 2 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
ifrunave = 1 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
nave = 5 ;five days run ave
ca = 40.0 ;day
cb = 8.0
rmv = True ;False;remove annaul cycle when bandpass filter

month_s = 6
month_e = 8
lagg = ispan(3,0,1) 
lagg = 3*lagg
nlag = dimsizes(lagg)

siglvl = 0.05
sig_option  = False;True ;if False, only the significant values are shaded 
uvsig_option= True ;False;if false, all wind will be drawed
vect_option = True ;False;
cont_option = False;True ;
cont_var  = (/"Z3" ,"intQ1"      ,"PSL","Div"                ,"OMEGA","U"  /)  ;this is the variable used to draw, can used as center string
cont_unit = (/"gpm","W/m~S~2~N~" ,"hPa","10~S~-6~N~s~S~-1~N~","Pa/s" ,"m/s"/)
nvc = 4
contsig  = True ;False;;if True, only the significant values are contoured 
barshare = True ;False;
barvert  = False;True ;

gtopfile= "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
years   = ispan(1979,2005,1)   ;the calculate time
path2   = "/home/ys17-19/data-observation"
case1   = (/"CN05.1","APHRO"/);read preci 
filein1 = path2 +(/"/CN05.1/CN05.1_Pre_1961_2018_daily_025x025.nc",\
                   "/APHRO_MA_025deg_V1003R1/APHRO_MA_025deg_V1003R1."/); + years(nt) + ".nc" ;1951-2007
case2   = (/"ERAI","NCEP1","NCEP2"/) ;read wind
filein2 = path2 +(/"/ERA-interim/pressure/",\ ;/uwnd/uwnd.interim.1993.nc 1979-2018
                   "/NCEP_NCAR/daily/",\ ;/uwnd/uwnd.1979.nc 1948-2018
                   "/NCEP_DOE/daily/"/)  ;/uwnd/uwnd.1979.nc,   short  1979-2020
nc1 = 0
nc2 = 0
pre_case   = (/"(a","(b","(c","(d","(a","(b","(c","(d","(i","(j","(k","(l","(i","(j","(k","(l","(i","(j","(k","(l"/)
fig_out    = "eps";"x11";"eps";
;fig_out    = "pdf";"x11";"eps";

fig_name = path + "/fig/" + case1(nc1) + "-project_ts_regr_" + lev(nl) + "wind_" + draw_var(nv)+"_"+cb+"-"+ca
if(cont_option) then 
fig_name = path + "/fig/" + case1(nc1) + "-project_ts_regr_" + lev(nl) + draw_var(nv) + "_" + cont_var(nvc)
end if

;area: east Asia
lats = (/0  ,0  ,0  /);-20
latn = (/90 ,80 ,70 /);60
lonl = (/0  ,0  ,40 /);20
lonr = (/360,180,160/);220
lat_spacing = 20
lon_spacing = 20
na = 2

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nyear = dimsizes(years)
;=====================================================
;read data and filter and regression
;==================================================
ft   = addfile("/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.PRECC.nc","r")
time = cd_calendar(ft->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.(month_s-1).and.time(:,1).le.(month_e+1))
end if

if(nc2.eq.0) then 
    f   = addfile(filein2(nc2)+"uwnd/uwnd.interim."+years(0)+".nc","r")
else
    f   = addfile(filein2(nc2)+"uwnd/uwnd."+years(0)+".nc","r")
end if
vars := f->uwnd(0,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
dim  := dimsizes(vars)
nlat := dim(0)
nlon := dim(1)
nvar = 3
ntime=dimsizes(select_time)
var  = new((/nvar,nyear,ntime,nlat,nlon/),float)

do nt = 0, nyear-1 ,1
if(isleapyear(years(nt))) then ;leap year,366 day
if(nc2.eq.0) then 
    inv = 0
    f   = addfile(filein2(nc2)+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->uwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc2)+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->vwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc2.eq.1) then 
    inv = 0
    f   = addfile(filein2(nc2)+"uwnd/uwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->uwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc2)+"vwnd/vwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->vwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc2.eq.2) then 
    inv = 0
    f   = addfile(filein2(nc2)+"uwnd/uwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = short2flt(f->uwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
    inv = 1
    f   = addfile(filein2(nc2)+"vwnd/vwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = short2flt(f->vwnd(select_time+1,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
end if
else
if(nc2.eq.0) then 
    inv = 0
    f   = addfile(filein2(nc2)+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc2)+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc2.eq.1) then 
    inv = 0
    f   = addfile(filein2(nc2)+"uwnd/uwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein2(nc2)+"vwnd/vwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(nc2.eq.2) then 
    inv = 0
    f   = addfile(filein2(nc2)+"uwnd/uwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = short2flt(f->uwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
    inv = 1
    f   = addfile(filein2(nc2)+"vwnd/vwnd."+years(nt)+".nc","r")
    var(inv,nt,:,:,:) = short2flt(f->vwnd(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}))
end if
end if
end do
if(nv.eq.0) then ;vort
    var(2,:,:,:,:) = uv2vr_cfd(var(0,:,:,:,:),var(1,:,:,:,:),var&lat,var&lon,2) 
    var(2,:,:,:,:) = var(2,:,:,:,:)*10^6
end if
if(nv.eq.1) then
    filein3 = "/home/ys17-19/data-observation/NOAA-OLR/daily/olr-daily_v01r02_"
    f = addfile(filein3+years(0)+"0101_"+years(0)+"1231.nc","r")
    vars1:= f->olr(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
    lat = vars1&lat
    lon = vars1&lon
    dim  = dimsizes(vars1)
    nlat1 = dim(0)
    nlon1 = dim(1)
    pvar := new((/nyear,ntime,nlat1,nlon1/),float)
    do nt = 0, nyear-1, 1
        f = addfile(filein3+years(nt)+"0101_"+years(nt)+"1231.nc","r")
        if(isleapyear(years(nt))) then ;leap year,366 day
        pvar(nt,:,:,:) = f->olr(select_time+1,{lats(na):latn(na)},{lonl(na):lonr(na)})
        else
        pvar(nt,:,:,:) = f->olr(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
        end if
    end do
    var(2,:,:,:,:) = linint2(lon,lat,pvar,False,var&lon,var&lat,0)
    ;the last number 0 is Reserved for future use.Currently not used, but set it to 0.
    delete([/vars1,pvar,lat,lon,nlat1,nlon1/])
end if
if(nv.eq.2) then 
    if(nc1.eq.0) then ;CN051
        f   = addfile(filein1(nc1),"r")
        vars1:= f->pre(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
        lat = vars1&lat
        lon = vars1&lon
    end if
    if(nc1.eq.1) then ;APHRO
        f   = addfile(filein1(nc1)+years(0)+".nc","r")
        vars1:= f->precip(0,{lats(na):latn(na)},{lonl(na):lonr(na)});mm/day
        lat = vars1&latitude
        lon = vars1&longitude
    end if
    dim  = dimsizes(vars1)
    nlat1 = dim(0)
    nlon1 = dim(1)
    pvar := new((/nyear,ntime,nlat1,nlon1/),float)
    if(nc1.eq.0) then 
        f   = addfile(filein1(nc1),"r")
        t   = cd_calendar(f->time,0)
        select_time0 = ind(t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1) \
                        .and.t(:,1).ge.(month_s-1).and.t(:,1).le.(month_e+1))
        pvar = reshape(f->pre(select_time0,{lats(na):latn(na)},{lonl(na):lonr(na)}),(/nyear,ntime,nlat1,nlon1/));mm/day
    end if
    if(nc1.eq.1) then 
        do nt = 0, nyear-1 ,1
        f  = addfile(filein1(nc1)+years(nt)+".nc","r")
        if(isleapyear(years(nt))) then ;leap year,366 day
        pvar(nt,:,:,:) = f->precip(select_time+1,{lats(na):latn(na)},{lonl(na):lonr(na)})
        else
        pvar(nt,:,:,:) = f->precip(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
        end if
        end do
    end if
    var(2,:,:,:,:) = linint2(lon,lat,pvar,False,var&lon,var&lat,0)
    ;the last number 0 is Reserved for future use.Currently not used, but set it to 0.
    delete([/vars1,pvar,lat,lon,nlat1,nlon1/])
end if

dimh = 2
if(ifrunave.eq.1) then 
opt1= -1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
var = runave_n_Wrap(var, nave, opt1, dimh)
var = var - conform(var,dim_avg_n(var,1),(/0,2,3,4/))  ;remove annual cycle
tsfile = path + "/data/"+case1(0)+"-eof_ts_"+month_s+month_e+"_"+var_name(2)+"_"+nave+"run.nc"
figtle = "runave="+nave+"day"
end if
if(ifrunave.eq.2) then 
fca = 1.0/ca
fcb = 1.0/cb
opt2= False
tsfile  = path + "/data/"+case1(0)+"-eof_ts_"+month_s+month_e+"_"+var_name(2)+"_"+cb+"-"+ca+".nc"
figtle = "bw_filter="+cb+"-"+ca+"day"
if(rmv) then 
var = var - conform(var,dim_avg_n(var,1),(/0,2,3,4/))  ;remove annual cycle
tsfile  = path + "/data/"+case1(0)+"-eof_ts_"+month_s+month_e+"_"+var_name(2)+"_"+cb+"-"+ca+"-rmv.nc"
figtle = "bw_filter="+cb+"-"+ca+"day rmv"
end if
var = bw_bandpass_filter(var,fca,fcb,opt2,dimh)
end if

if(ifrunave2.eq.2) then 
    fca = 1.0/ca
    fcb = 1.0/cb
    opt2= False
    var = bw_bandpass_filter(var,fca,fcb,opt2,dimh)
end if

time := cd_calendar(ft->time(select_time),0)
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
;var := reshape(var(:,:,select_time,:,:),(/nvar,nyear*ntime,nlat,nlon/))
printVarSummary(var)

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

;define the thickness and height of the title,tick and label
if(na.eq.0) then 
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.02
   vcref_font     = 0.008
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.5
   font_index     = 22
else
   tick_font      = 0.025
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.025
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.7 ;3.0;
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
    res@mpGeophysicalLineColor      = "lightsalmon3";"black";"white"

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
    resvc@vcPositionMode            = "ArrowTail"  ;格点位置对应箭头的尾部
    resvc@vcGlyphStyle              = "CurlyVector";"LineArrow"  ;"Fillarrow"
    ;resvc@vcMinMagnitudeF           = 0.3          ; 小于该数值则不绘制
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = mp_thick-0.4 
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
if(na.eq.0) then 
    resvc@vcMinDistanceF            = 0.009
    resvc@vcRefLengthF              = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
else
    resvc@vcMinDistanceF            = 0.023
    resvc@vcRefLengthF              = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.008;0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.008;0.015
end if
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    resvc@vcRefAnnoOn                = False;True;
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.8          ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.4            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False;True
    resvc@vcRefAnnoString2On         = True ;False;
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor          = "grey"

;----------Set res line for contour-----------
    res_line                  = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False       ; Fill contours
    res_line@cnLinesOn        = True      ; turn off contour lines
    res_line@cnLineLabelsOn   = False      ; turn off contour line labels
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 5.0
    res_line@cnMonoLineColor  = True ;False
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 2
    res_line@gsnContourPosLineDashPattern = 0
    res_line@cnLevelSelectionMode  = "AutomaticLevels"   ; set manual contour levels
    ;res_line@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res_line@cnLevelSpacingF       = 1.0 ;        ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = mp_thick + 0.5 
    res_gtop@cnLineColor      = "lightsalmon3";"black"
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
    resp@gsnPanelDebug = True
    ;resp@gsnPanelScalePlotIndex     = 2 
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = True
    resp@lbTitleOn    =  True                ; turn on title
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.7
    else
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.95 ;0.05
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
    resp@lbTitleString            = unit(nv) 
    resp@lbTitlePosition          = "Right"              ; title position
    resp@lbTitleFontHeightF       = .015                ; make title smaller
    resp@lbTitleFont              = font_index
    resp@lbTitleOffsetF           = 0.01
    resp@lbTitleDirection         = "Across"             ; title direction

if(na.eq.0) then ; 850,700,600,500,200,100,10
    vcref         = (/5  ,5 ,5 ,5 /)  ;the level of 850hpa    ctrl,nudg,diff
else                               
    vcref         = (/0.5, 3, 1,3 /)  ;the level of 850hpa    ctrl,nudg,diff
end if

minlvl = (/-4.0,- 8 ,-2.4/)  ;q1 ;-80,-80
maxlvl = (/ 4.0,  8 , 2.4/)      ; 80, 80
spacig = (/ 0.5,  1 , 0.3/)      ; 10, 10

res@cnFillPalette = "BlueDarkRed18"
res@cnFillColors  = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/)

tsf = addfile(tsfile,"r")
eval = tsf->eval(nc1,:)
pcvar= tsf->pcvar(nc1,:)
np = -1
do ne = 0, 3, 1
ts = tsf->eof_ts(nc1,ne,:)
if(ifrunave2.eq.2) then 
    fca = 1.0/ca
    fcb = 1.0/cb
    opt2= False
    ts = bw_bandpass_filter(ts,fca,fcb,opt2,0)
end if

printVarSummary(ts)
regr := new((/nvar,nlag,nlat,nlon/),float)
tval := new((/nvar,nlag,nlat,nlon/),float)
df   := new((/nvar,nlag,nlat,nlon/),float)
prob := new((/nvar,nlag,nlat,nlon/),float)
do nlg = 0, nlag-1, 1
    rvar := reshape(var(:,:,select_time-lagg(nlg),:,:),(/nvar,nyear*ntime,nlat,nlon/))
    regr(:,nlg,:,:) = regCoef_n(ts,rvar,0,1)
    tval(:,nlg,:,:) = onedtond(regr@tval ,(/nvar,nlat,nlon/))
    df(:,nlg,:,:)   = onedtond(regr@nptxy,(/nvar,nlat,nlon/))
end do
if(ifrunave.eq.1) then 
df = toint(df/nave)
end if
copy_VarCoords(vars, regr(0,0,:,:))
b := tval 
b = 0.5
prob  := betainc( df/(df+tval^2), df/2.0, b )
delete([/tval,df/])

if(uvsig_option) then 
notsigu = mask(regr(0,:,:,:),(prob(0,:,:,:).lt.siglvl).or.(prob(1,:,:,:).lt.siglvl),False)
notsigv = mask(regr(1,:,:,:),(prob(0,:,:,:).lt.siglvl).or.(prob(1,:,:,:).lt.siglvl),False)
regr(0,:,:,:) = mask(regr(0,:,:,:),(prob(0,:,:,:).lt.siglvl).or.(prob(1,:,:,:).lt.siglvl),True)
regr(1,:,:,:) = mask(regr(1,:,:,:),(prob(0,:,:,:).lt.siglvl).or.(prob(1,:,:,:).lt.siglvl),True)
copy_VarCoords(vars, notsigu(0,:,:))
copy_VarCoords(vars, notsigv(0,:,:))
printVarSummary(notsigu)
end if
if(.not.(sig_option)) then
regr(2,:,:,:) = mask(regr(2,:,:,:),prob(2,:,:,:).lt.siglvl,True)
end if
prob   = mask(prob,prob.lt.siglvl,True)
printVarSummary(regr)

nne = ne+1
do nlg = 0,nlag-1,1
    ;np = nlg*4+ne
    np = np + 1
    if(lagg(nlg).eq.0) then 
    res@gsnLeftString   = pre_case(np) + ") OBS (" + lagg(nlg) + " day)" ; 
    else                                                    
    res@gsnLeftString   = pre_case(np) + ") OBS (-" + lagg(nlg) + " day)";  
    end if
    ;if(nlg.eq.0) then 
    ;res@tmYLLabelsOn         = True 
    ;else
    ;res@tmYLLabelsOn         = False
    ;end if
    ;if(ne.eq.3)
    ;res@tmXBLabelsOn         = True 
    ;else
    ;res@tmXBLabelsOn         = False
    ;end if
    res@gsnRightString  ="" 
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np)    = gsn_csm_contour_map_ce(wks,regr(2,nlg,:,:),res)
    if(sig_option) then 
    sig(np) = gsn_csm_contour(wks,prob(2,nlg,:,:),res_sig)
    overlay(plot(np),sig(np))
    end if

    if(vect_option) then
    if(nlg.eq.nlag-1.and.ne.eq.3) then
    resvc@vcRefAnnoOn                = True ;False;
    else
    resvc@vcRefAnnoOn                = False;True;
    end if
    resvc@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl)+"m/s"
    resvc2@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc2@vcRefAnnoString2     = vcref(nl)+"m/s"
    if(uvsig_option) then
    wind2(np)   = gsn_csm_vector(wks,notsigu(nlg,:,:),notsigv(nlg,:,:),resvc2)
    overlay(plot(np),wind2(np))
    end if
    wind(np)    = gsn_csm_vector(wks,regr(0,nlg,:,:),regr(1,nlg,:,:),resvc)
    overlay(plot(np),wind(np))
    end if
    
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    if(cont_option)
    res_line@cnLevelSpacingF    = cont_spacing(nc) 
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    plot2(np)   = gsn_csm_contour(wks,regr(3,nlg,:,:),res_line)
    overlay(plot(np),plot2(np))
    end if
    
    print(max(regr(2,nlg,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    print(min(regr(2,nlg,{lats(na):latn(na)},{lonl(na):lonr(na)})))
end do 
end do
resp@txString = case2(nc2)+" "+figtle+" "+month_s+"-"+month_e+" "+lev(nl) + "hPa " + draw_var(nv) + "(" + unit(nv) + ") " ;
gsn_panel(wks,plot,(/4,4/),resp)
end

