;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
var_option = 0  ;0 means Qd; 1 means f_Qd; 2 means dzdt; 3 means duTdx; 4 means dzdt induced by b1
area_option = 1

if(var_option.eq.0) then 
    draw_var = (/"sigm1"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"sigm1"/)  ;
    unit     = (/"10~S~-4~N~J/(Pa~S~2~N~Kg)"/)  ;this can be used as right string 
    nv = 0
end if
;months   = (/"DJF"/)
;lev  = (/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100/) 
months   = (/"JF","MA","MJ","JA"/)
;lev  = (/1000, 925, 850, 700, 600, 500,400,300,250,200,175,150,125,100/) 
;lev  = (/1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100/)
;lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/)
lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100/)
;lev  = (/800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100/)

if(area_option.eq.0) then 
area = (/"TP","SEA","NEA","NH","CA"/)
lats = (/15  , 20  , 30  , 0  ,10  /);20 ;-20
latn = (/55  , 30  , 45  , 89 ,20  /) ;55 ;60
lonl = (/70  ,110  ,110  , 1  ,0   /) ;30 ;20
lonr = (/105 ,130  ,130  ,359 ,50  /) ;150;220
end if
if(area_option.eq.1) then 
area = (/"NH","NH","NH","NH","NH"/)
lats = (/ 0  , 20 , 30 , 50 , 60 /);20 ;-20
latn = (/ 89 , 89 , 89 , 89 , 89 /) ;55 ;60
lonl = (/ 1  , 1  , 1  , 1  , 1  /) ;30 ;20
lonr = (/359 ,359 ,359 ,359 ,359 /) ;150;220
end if
month_s = 1
month_e = 8

;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
case   = "ERA-Interim"
;case   = "NCEP1"

if(var_option.le.0)
    filename = path + "mdata/" + case + "_Clim_month_Qd_Qeddy_A_sigm.nc"
    fig_name = path + "fig/1p2x3_height_lat_" + var_name(0)
end if

;fig_title  = "Seasonal " + draw_var(nv) + " variance in " ;+season
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
fig_out    ="pdf";"x11"

narea  = dimsizes(area)
nmonth = dimsizes(months)
nlev   = dimsizes(lev)
var_m  = new((/narea,nmonth,nlev/),double)
;var_m := new((/nvar,nmonth,nlev,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/nmonth/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.015
   label_font     = 0.015
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
pi = atan(1.0)*4
f = addfile(filename,"r")
do na = 0, narea-2,1
    var := f->$var_name(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    var_m(na,:,:) = wgt_areaave(var, cos(var&lat*pi/180.0), 1.0, 0) ;If opt = 0, the area average is calculated using available non-missing data
end do
copy_VarMeta(var(0,:,0,0),var_m(0,0,:))
printVarSummary(var_m)

if(var_option.eq.0) then 
    var_m = var_m*10^4
end if

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

    resxy     = True
    resxy@vpWidthF  = 0.4
    resxy@vpHeightF = 0.3
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
   
    if(area_option.eq.0) then 
    resxy@xyDashPatterns    = (/0,0,0,0,0/) ;,0,0
    resxy@xyLineThicknesses = (/4.0, 4.0, 4.0, 4.0, 4.0/) ;, 4.0, 4.0  ; make second line thicker
    resxy@xyLineColors      = (/"black","red","blue","grey","orange"/); ,"darkorange3","chartreuse4"
    resxy@xyExplicitLegendLabels = area 
    end if
    if(area_option.eq.1) then 
    resxy@xyDashPatterns    = (/0,0,0,0,0/) ;,0,0
    resxy@xyLineThicknesses = (/4.0, 4.0, 4.0, 4.0, 4.0/) ;, 4.0, 4.0  ; make second line thicker
    resxy@xyLineColors      = (/"black","red","blue","grey","orange"/); ,"darkorange3","chartreuse4"
    ;resxy@xyDashPatterns    = (/0,0,0/) ;,0,0
    ;resxy@xyLineThicknesses = (/4.0, 4.0, 4.0/) ;, 4.0, 4.0  ; make second line thicker
    ;resxy@xyLineColors      = (/"black","red","blue"/); ,"darkorange3","chartreuse4"
    resxy@xyExplicitLegendLabels = lats + " - " + latn + "N"
    end if
    resxy@tiXAxisOn = False 
    resxy@tiXAxisString = "lag(days)"
    resxy@tiYAxisOn = False
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font
    resxy@tmXBLabelFontThicknessF  = tick_thick 
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    ;resxy@tmXBMode           = "Manual"
    resxy@tmXBTickSpacingF   = 0.5
    resxy@tmXBTickStartF     = 0
    resxy@tmXBTickEndF       = 1.5
    resxy@trXMinF = 0 ;floor(min(var)); -0.4;
    resxy@trXMaxF = 2 ;ceil(max(var)) ; 0.2 ;
    ;resxy@tmXBMode            = "Explicit"
    ;resxy@tmXBValues          = fspan(0,1.5,5)
    ;resxy@tmXBLabels          = fspan(0,1.5,5)
    
    resxy@trYReverse        = True
    resxy@xyYStyle          = "Log"
    resxy@xyComputeYMin     = True
    resxy@tmYLMode          = "Explicit"                ; explicit labels
    resxy@tmYLValues        = (/ 1000, 700,500,400,300,200,100/)
    resxy@tmYLLabels        = ""+resxy@tmYLValues         ; make strings
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    ;resxy@gsnYRefLine           = 0.0
    ;resxy@gsnYRefLineColor      = "black"
    ;resxy@gsnYRefLineThicknessF = 5.0
    ;resxy@gsnYRefLineDashPattern= 0

    ;legend
    resxy@pmLegendDisplayMode = "NoCreate";"Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    ;resxy@xyExplicitLegendLabels = label
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.09
    resxy@pmLegendHeightF  = 0.3
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.05   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.3    ; move ref vector along X-axis
    
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF   = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = "" 
    resxy@gsnRightString   = ""

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do nm = 0,nmonth-1,1
    np = nm
    if(nm.eq.3) then 
    resxy@pmLegendDisplayMode = "Always"
    else
    resxy@pmLegendDisplayMode = "NoCreate";"Always"
    end if
    resxy@gsnLeftString   = pre_case(np) + ") " + case + " " + months(nm) 
    resxy@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")"
    plot(np) = gsn_csm_xy(wks,var_m(:,nm,:),var_m&lev,resxy)
    print(max(var_m(:,nm,:)))
    print(min(var_m(:,nm,:)))
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/2,2/),resp)
end

