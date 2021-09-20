;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
path    = "/home/ys17-19/renql/project/2021variability"
var_name = (/"vort","preci"/)
unit     = (/"10~S~-6~N~ m~S~-1~N~","mm/day"/)
draw_var = (/"vort","preci"/)      
lev = (/850,500,200,700,600,  100,10/)
nl = 2
nv = 1

ifrunave2= 2 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
ifrunave = 1 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
nave = 5 ;five days run ave
ca = 40.0 ;day
cb = 8.0
rmv = False;True ;remove annaul cycle when bandpass filter

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

case1   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
case2   = (/"CTL" ,"NT"   ,"NM"    ,"NMT" /) 
gtopfile= "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
years   = ispan(1979,2005,1)   ;the calculate time
path2   = "/home/ys17-19/renql"
filein  = path2 +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                   "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                   "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                   "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
pre_case   = (/"(b","(c","(d","(e","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t"/)
fig_out    = "eps";"x11";"eps";
;fig_out    = "pdf";"x11";"eps";

fig_name = path + "/fig/" + case1(0) + "-project_ts_regr_" + lev(nl) + "wind_" + draw_var(nv)+"_"+cb+"-"+ca
if(cont_option) then 
fig_name = path + "/fig/" + case1(0) + "-project_ts_regr_" + lev(nl) + draw_var(nv) + "_" + cont_var(nvc)
end if

;area: east Asia
lats = (/0  ,0  /);-20
latn = (/90 ,70 /);60
lonl = (/0  ,40  /);20
lonr = (/360,160/);220
lat_spacing = 20
lon_spacing = 30
na = 1

f    = addfile(filein(0)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.(month_s-1).and.time(:,1).le.(month_e+1))
end if

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ntime = dimsizes(select_time)
ncase = dimsizes(case1)
nyear = dimsizes(years)
nlat  = dimsizes(phis&lat)
nlon  = dimsizes(phis&lon)
nvar  = 3
var   = new((/ncase,nvar,nyear,ntime,nlat,nlon/),float)

;=====================================================
;read data and filter and regression
;==================================================
do nc = 0, ncase-1, 1
do nt = 0, nyear-1 ,1
    inv = 0
    f   = addfile(filein(nc)+years(nt)+".daily.U.nc","r")
    var(nc,inv,nt,:,:,:) = f->U(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    f   = addfile(filein(nc)+years(nt)+".daily.V.nc","r")
    var(nc,inv,nt,:,:,:) = f->V(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    if(nv.eq.1) then 
    inv = 2
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    var(nc,inv,nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    var(nc,inv,nt,:,:,:) = var(nc,inv,nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    end if
    if(nv.ge.2) then 
    inv = 2
    f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var(nc,inv,nt,:,:,:) = f->$var_name$(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    end if
end do
end do
if(nv.eq.0) then ;vort
    var(:,2,:,:,:,:) = uv2vr_cfd(var(:,0,:,:,:,:),var(:,1,:,:,:,:),var&lat,var&lon,2) 
    var(:,2,:,:,:,:) = var(:,2,:,:,:,:)*10^6
end if
if(nv.eq.1) then ;preci
    var(:,2,:,:,:,:) = var(:,2,:,:,:,:)*1000*3600*24
end if

dimh = 3
if(ifrunave.eq.1) then 
opt1= -1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
var = runave_n_Wrap(var, nave, opt1, dimh)
var = var - conform(var,dim_avg_n(var,2),(/0,1,3,4,5/))  ;remove annual cycle
tsfile  = path + "/data/"+case1(0)+"-eof_project_ts_"+month_s+month_e+"_"+var_name(nv)+"_"+nave+"run.nc"
figtle = "runave="+nave+"day"
end if
if(ifrunave.eq.2) then 
fca = 1.0/ca
fcb = 1.0/cb
opt2= False
tsfile  = path + "/data/"+case1(0)+"-eof_project_ts_"+month_s+month_e+"_"+var_name(nv)+"_"+cb+"-"+ca+".nc"
figtle = "bw_filter="+cb+"-"+ca+"day"
if(rmv) then 
var = var - conform(var,dim_avg_n(var,2),(/0,1,3,4,5/))  ;remove annual cycle
tsfile  = path + "/data/"+case1(0)+"-eof_project_ts_"+month_s+month_e+"_"+var_name(nv)+"_"+cb+"-"+ca+"-rmv.nc"
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

time := cd_calendar(f->time(select_time),0)
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
;var := reshape(var(:,:,:,select_time,:,:),(/ncase,nvar,nyear*ntime,nlat,nlon/))

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
   mp_thick       = 1.7;3.0;
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
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = mp_thick-0.4 
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "grey20"
if(na.eq.0) then 
    resvc@vcMinDistanceF            = 0.009
    resvc@vcRefLengthF              = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
else
    resvc@vcMinDistanceF            = 0.025
    resvc@vcRefLengthF              = 0.05 ;0.05           ; define length of vec ref
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
    resp@pmLabelBarWidthF  = 0.9 ;0.05
    resp@pmLabelBarHeightF = 0.04;0.8
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
    resp@lbLabelFontHeightF = 0.012 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    resp@lbTitleString            = unit(nv) 
    resp@lbTitlePosition          = "Right"              ; title position
    resp@lbTitleFontHeightF       = .012                 ; make title smaller
    resp@lbTitleFont              = font_index
    resp@lbTitleOffsetF           = 0.012
    resp@lbTitleDirection         = "Across"             ; title direction

if(na.eq.0) then ; 850,700,600,500,200,100,10
    vcref         = (/5  ,5 ,5 ,5 /)  ;the level of 850hpa    ctrl,nudg,diff
else                               
    vcref         = (/0.5, 3, 1,3 /)  ;the level of 850hpa    ctrl,nudg,diff
end if

minlvl = (/-4.0,-3.5/)  ;q1 ;-80,-80
maxlvl = (/ 4.0, 3.5/)      ; 80, 80
spacig = (/ 0.5, 0.5/)      ; 10, 10

res@cnFillPalette = "BlueDarkRed18"
res@cnFillColors  = (/0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17/)

tsf = addfile(tsfile,"r")
eval = tsf->eval
pcvar= tsf->pcvar
pcvar= round(pcvar*10,1)/10.0 ;Round to one decimal place
do ne = 0,0, 1;5, 1
ts = tsf->eof_ts(2:5,ne,:)
if(ifrunave2.eq.2) then 
    fca = 1.0/ca
    fcb = 1.0/cb
    opt2= False
    ts = bw_bandpass_filter(ts,fca,fcb,opt2,1)
end if

printVarSummary(ts)
regr := new((/ncase,nvar,nlag,nlat,nlon/),float)
tval := new((/ncase,nvar,nlag,nlat,nlon/),float)
df   := new((/ncase,nvar,nlag,nlat,nlon/),float)
prob := new((/ncase,nvar,nlag,nlat,nlon/),float)
do nc = 0, ncase-1, 1
do nlg = 0, nlag-1, 1
    rvar := reshape(var(:,:,:,select_time-lagg(nlg),:,:),(/ncase,nvar,nyear*ntime,nlat,nlon/))
    regr(nc,:,nlg,:,:) = regCoef_n(ts(nc,:),rvar(nc,:,:,:,:),0,1)
    tval(nc,:,nlg,:,:) = onedtond(regr@tval ,(/nvar,nlat,nlon/))
    df(nc,:,nlg,:,:)   = onedtond(regr@nptxy,(/nvar,nlat,nlon/))
end do
end do
copy_VarCoords(phis, regr(0,0,0,:,:))
if(ifrunave.eq.1) then 
df = toint(df/nave)
end if
b := tval 
b = 0.5
prob  := betainc( df/(df+tval^2), df/2.0, b )
delete([/tval,df/])

if(uvsig_option) then 
notsigu = mask(regr(:,0,:,:,:),(prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl),False)
notsigv = mask(regr(:,1,:,:,:),(prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl),False)
regr(:,0,:,:,:) = mask(regr(:,0,:,:,:),(prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl),True)
regr(:,1,:,:,:) = mask(regr(:,1,:,:,:),(prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl),True)
copy_VarCoords(phis, notsigu(0,0,:,:))
copy_VarCoords(phis, notsigv(0,0,:,:))
printVarSummary(notsigu)
end if
if(.not.(sig_option)) then
regr(:,2,:,:,:) = mask(regr(:,2,:,:,:),prob(:,2,:,:,:).lt.siglvl,True)
end if
prob   = mask(prob,prob.lt.siglvl,True)
printVarSummary(regr)

np = -1
nne = ne+1
resp@txString = figtle+" month: "+month_s+"-"+month_e+" EOF"+nne+": "+eval(ne)+"  "+pcvar(ne)+"%" 
do nc = 0,ncase-1,1
do nlg = 0,nlag-1,1
    np = np + 1
    if(lagg(nlg).eq.0) then 
    res@gsnLeftString   = case2(nc) +" (" + lagg(nlg) + " day)";  
    else                                                                               
    res@gsnLeftString   = case2(nc) +" (-" + lagg(nlg) + " day)";  
    end if
    res@gsnRightString  = "" 
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np)    = gsn_csm_contour_map_ce(wks,regr(nc,2,nlg,:,:),res)
    if(sig_option) then 
    sig(np) = gsn_csm_contour(wks,prob(nc,2,nlg,:,:),res_sig)
    overlay(plot(np),sig(np))
    end if

    if(vect_option) then
    if(nlg.eq.nlag-1.and.nc.eq.ncase-1) then
    resvc@vcRefAnnoOn                = True ;False;
    else
    resvc@vcRefAnnoOn                = False;True;
    end if
    resvc@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl)+"m/s"
    resvc2@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc2@vcRefAnnoString2     = vcref(nl)+"m/s"
    ;if(uvsig_option) then
    ;wind2(np)   = gsn_csm_vector(wks,notsigu(nc,nlg,:,:),notsigv(nc,nlg,:,:),resvc2)
    ;overlay(plot(np),wind2(np))
    ;end if
    wind(np)    = gsn_csm_vector(wks,regr(nc,0,nlg,:,:),regr(nc,1,nlg,:,:),resvc)
    overlay(plot(np),wind(np))
    end if
    
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    if(cont_option)
    res_line@cnLevelSpacingF    = cont_spacing(nc) 
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    plot2(np)   = gsn_csm_contour(wks,regr(nc,3,nlg,:,:),res_line)
    overlay(plot(np),plot2(np))
    end if
    
    if(nc.eq.3) then
    line(np)  = gsn_add_polyline( wks, plot(np), lonx(:,0), latx(:,0), lnres )
    line1(np) = gsn_add_polyline( wks, plot(np), lonx(:,1), latx(:,1), lnres )
    end if
    if(nc.eq.1.or.nc.eq.2) then
    line(np) = gsn_add_polyline( wks, plot(np), lonx(:,nc-1), latx(:,nc-1), lnres )
    end if
    
    print(max(regr(nc,2,nlg,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    print(min(regr(nc,2,nlg,{lats(na):latn(na)},{lonl(na):lonr(na)})))
end do 
end do 
gsn_panel(wks,plot,(/4,4/),resp)
end do
end

