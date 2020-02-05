;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
path    = "/home/ys17-19/renql/project/TP_NUDG"
case_option = 1 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h","NUDG24h"/)
    case2   = (/"CTRL","NG6h-CTRL","NG24h-CTRL"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
end if
if(case_option.eq.1) then
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
end if
if(case_option.eq.2) then
    case1   = (/"F2000","NUDG"/)
    case2   = (/"F2000","NUDG-F2000"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 30
end if

;Q1 is calculated from temperature equation; q1 is calculated from (QRL,QRS,DTCOND,DTV)
draw_var = (/"intQ1","preci" ,"OMEGA"          ,"Z3" ,"intq1"/)  ;this is the variable used to draw, can used as center string
unit     = (/"W/kg" ,"mm/day","10~S~-2~N~ Pa/s","gpm","W/kg"/)    ;this can be used as right string 
draw_varname = (/"Q1","preci","omega","Z3","q1"/)  ;this is the variable used to draw, can used as center string
lev = (/850,700,500,200/)
nl = 1
nv = 3
diff = False;True

var_option = 0;2;1
if(var_option.eq.0) then 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
imonth = ispan(0,3,1)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-" + lev(nl) + "wind_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
end if
if(var_option.eq.1) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;imonth = ispan(0,3,1);(/4,5,6,7/)
;imonth = ispan(4,7,1);(/4,5,6,7/)
imonth = ispan(8,11,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-" + lev(nl) + "wind_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
end if
if(var_option.eq.2) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"JF","MA","MJ","JA"/)
imonth = ispan(2,3,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-2month_" + lev(nl) + "wind_" + draw_var(nv)
end if

;pre_case   = (/"(e","(f","(g","(h","(e","(f","(g","(h","(i","(j","(k","(l"/)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    = "pdf";"x11""eps";

;area: east Asia
lats = (/0  ,0  ,0  ,-70,-70,-70,10 ,-70, 0/);-20
latn = (/89 ,70 ,70 , 0 , 0 , 0 ,55 ,70 ,50 /);60
lonl = (/0  ,150,270,30 ,150,270,90 ,0  ,40 /);20
lonr = (/360,270,390,150,270,390,150,360,130/);220
siglvl = 0.05
lat_spacing = 15
lon_spacing = 20
na = 8

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case1)
nmonth = dimsizes(imonth)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nvar = 3
var    = new((/ncase  ,nvar,nmonth,nlat,nlon/),float)
vari   = new((/ncase  ,nvar,nmonth,nlat,nlon/),float)
prob   = new((/ncase-1,nvar,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot   = new((/ncase*nmonth/),graphic)
    wind   = new((/ncase*nmonth/),graphic)
    wind2  = new((/ncase*nmonth/),graphic)
    gtop   = new((/ncase*nmonth/),graphic)
    line   = new((/ncase*nmonth/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
if(na.eq.0) then 
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   vcref_font     = 0.01
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 1.5
   font_index     = 22
else
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.0
   mp_thick       = 2.5
   font_index     = 22
end if

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-1, 1
inv = 0
f = addfile(filename(nf)+"U.nc","r")
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
inv = 1
f = addfile(filename(nf)+"V.nc","r")
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
inv = 2
f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
if(nv.le.1.or.nv.eq.4) then
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
else
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if
end do

if(diff) then 
vari = where(abs(vari).le.0.000001,0.1,vari)
do nf = 0 ,ncase-2,1
prob(nf,:,:,:,:) = ttest(var(0,:,:,:,:),vari(0,:,:,:,:),nsample,var(nf+1,:,:,:,:),vari(nf+1,:,:,:,:),nsample,True,False)
var(nf+1,:,:,:,:)= var(nf+1,:,:,:,:) - var(0,:,:,:,:)   ;get the difference between CTRL and NUDG
end do
end if

if(nv.eq.2) then
    var(:,2,:,:,:) = var(:,2,:,:,:)*100 ;omega*100
end if

if(nv.eq.4) then
    cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
    var(:,2,:,:,:) = var(:,2,:,:,:)*cp ;Q1,transform unit from K/S to W/kg
end if

if(diff) then 
uvar = mask(var(1:(ncase-1),0,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),False)
vvar = mask(var(1:(ncase-1),1,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),False)
var(1:(ncase-1),0,:,:,:) = mask(var(1:(ncase-1),0,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),True)
var(1:(ncase-1),1,:,:,:) = mask(var(1:(ncase-1),1,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),True)
;var(1:(ncase-1),2,:,:,:) = mask(var(1:(ncase-1),2,:,:,:),prob(:,2,:,:,:).lt.siglvl,True)
;prob = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(var(0,0,0,:,:),prob(0,0,0,:,:))
copy_VarMeta(var(0,0,0,:,:),uvar(0,0,:,:))
copy_VarMeta(var(0,0,0,:,:),vvar(0,0,:,:))
printVarSummary(prob)
end if
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

;do na = 0, 5, 1
print("draw lats = " + lats(na) + ", latn = " + latn(na) + ", lonl = " + lonl(na) + ", lonr = " + lonr(na))
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "landscape" ;"portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
if(na.eq.2.or.na.eq.5) then
    res@mpCenterLonF     = 0.
else
    res@mpCenterLonF     = 180.
end if
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
    ;res@cnLineThicknessF = 0
    ;res@gsnContourZeroLineThicknessF = 4.0
    ;res@cnLineColor      = "blue"
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 1          ; set min contour level
    ;res@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    ;res@pmLabelBarWidthF  = 0.65
    ;res@pmLabelBarHeightF = 0.08

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
if(na.eq.0) then 
    resvc@vcMinDistanceF             = 0.01
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
else
    resvc@vcMinDistanceF             = 0.025
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.015
end if
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.3           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor           = "grey"

;----------Set res line-----------
    res_line                  = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False       ; Fill contours
    res_line@cnLinesOn        = True      ; turn off contour lines
    res_line@cnLineLabelsOn   = False      ; turn off contour line labels
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 2.5
    res_line@cnMonoLineColor  = True ;False
    res_line@cnLineColor      = "black"
    res_line@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_line@cnLevels         = (/20/)
    res_line@lbLabelBarOn         = False   ; the label bar
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 2.0
    ;lnres@gsLineLabelString    = "TP"

    lnres2                      = True
    lnres2@gsLineColor          = "black"
    lnres2@gsLineDashPattern    = 1
    lnres2@gsLineThicknessF     = 3.5
    ;lnres@gsLineLabelString    = "Westly"

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = mp_thick
    res_gtop@cnLineColor      = "black"
    res_gtop@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_gtop@cnLevels         = (/1500/)
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
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = False; True
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.05
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(na.eq.0) then 
    vcref         = (/(/10,5 ,5 /), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/10,5 ,5 /), \  ;the level of 700hpa    ctrl,nudg,diff
                      (/7 ,3 ,3 /), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/20,10,10/)/)   ;the level of 200hpa    ctrl,nudg,diff
else
    vcref         = (/(/7 ,7 ,2 ,2 /), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/7 ,7 ,2 ,2 /), \  ;the level of 700hpa    ctrl,nudg,diff
                      (/7 ,7 ,3 ,3 /), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/15,15,4 ,4 /)/)   ;the level of 200hpa    ctrl,nudg,diff
end if

if(nv.eq.0) 
    minlvl = (/0  ,-160,-160/)  ;q1 ;-80,-80
    maxlvl = (/300, 160, 160/)      ; 80, 80
    spacig = (/20 ,  20,  20/)      ; 10, 10
end if
if(nv.eq.4) 
    minlvl = (/0  ,-160,-160/)  ;q1 ;-80,-80
    maxlvl = (/300, 160, 160/)      ; 80, 80
    spacig = (/20 ,  20,  20/)      ; 10, 10
end if

if(nv.eq.1) then ;preci
    minlvl = (/ 1 ,-4.0,-4.0/);-8,-8,-8,-8
    maxlvl = (/16 ,4.0 ,4.0 /);8 ,8 ,8 ,8 
    spacig = (/ 1 ,0.5 ,0.5 /);1 ,1 ,1 ,1 
end if

if(nv.eq.2) 
    minlvl = (/-8 ,-4 , -4/)  ;omega*100
    maxlvl = (/ 8 , 4 ,  4/)
    spacig = (/ 1 ,0.5,0.5/)
end if

if(nv.eq.3) ;Z3
if(nl.eq.0) then 
    minlvl = (/1300,1300,-40,-40/) 
    maxlvl = (/1620,1620, 40, 40/)
    spacig = (/ 20 , 20 , 5 , 5 /) 
end if
if(nl.eq.1) then 
    minlvl = (/2660,2660,-40,-40/) 
    maxlvl = (/3140,3140, 40, 40/)
    spacig = (/ 30 , 30 , 5 , 5 /) 
end if
if(nl.eq.2) then 
    minlvl = (/5730,-16,-16/)
    maxlvl = (/5880, 16, 16/)
    spacig = (/  10, 2 , 2 /) 
end if
if(nl.eq.3) then 
    minlvl = (/12063,-40,-40/) 
    maxlvl = (/12558, 40, 40/)
    spacig = (/  33 , 5 , 5 /) 
end if
end if

do nm = 0,nmonth-1,1
do nc = 0, ncase-1,1
if(nc.gt.2.or.nv.eq.2) then
    res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
else
   res@cnFillPalette = "precip3_16lev"
   ;blue and yellow and red, has 17 colors,2-18,2 is white
end if 
    ;do nm = 0,nmonth-1,1
        np = nc
        ;np = nmonth*nc + nm
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc)
        res@gsnRightString  = months(imonth(nm)) + " " +  draw_varname(nv) + " (" + unit(nv) + ") " + lev(nl) + "hPa"
        ;res@gsnRightString  = draw_varname(nv) + " (" + unit(nv) + ")" 
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        resvc@vcRefMagnitudeF      = vcref(nl,nc)             ; define vector ref mag
        resvc@vcRefAnnoString2     = vcref(nl,nc)+"m/s"
        resvc2@vcRefMagnitudeF      = vcref(nl,nc)             ; define vector ref mag
        resvc2@vcRefAnnoString2     = vcref(nl,nc)+"m/s"
        plot(np)    = gsn_csm_contour_map_ce(wks,var(nc,2,nm,:,:),res)
        if(diff.and.nc.gt.0)
        wind2(np)   = gsn_csm_vector(wks,uvar(nc-1,nm,:,:),vvar(nc-1,nm,:,:),resvc2)
        overlay(plot(np),wind2(np))
        end if
        wind(np)    = gsn_csm_vector(wks,var(nc,0,nm,:,:),var(nc,1,nm,:,:),resvc)
        overlay(plot(np),wind(np))
        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        if(na.eq.0.or.na.eq.6) then
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        end if
        print(max(var(nc,2,nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        print(min(var(nc,2,nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        ;draw(plot(np))
        ;frame(wks)
    end do 
    gsn_panel(wks,plot,(/ncase,1/),resp)
;end do ;for the month 

    ;res1 = True
    ;res2 = True
    ;res1@gsnMaximize = True
    ;attach1 = gsn_attach_plots(plot(0),plot(1:3),res1,res2)
    ;attach2 = gsn_attach_plots(plot(4),plot(5:7),res1,res2)
    ;attach3 = gsn_attach_plots(plot(8),plot(9:11),res1,res2)
    ;gsn_panel(wks,plot(0:8:4),(/3,1/),resp)
;    resp@txString      = fig_title
    ;gsn_panel(wks,plot,(/2,2/),resp)
    ;gsn_panel(wks,plot,(/4,1/),resp)
end do
;    gsn_panel(wks,plot,(/ncase,nmonth/),resp)
end

