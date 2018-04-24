;***************************************************************************
;filename:      draw_1p3x2_preci_JJAts.ncl
;plot the JJA times of every year's annual cycle(id=151~242) of TP average UVT(27 curves) 
;then overlaying the climate annual cycle(black bold curve)
;and one panel has six plots.
;                                            by Ql_Ren
;                                           2017/03/29
;******************************************************************************
begin
draw_var = (/"preci"/)  ;this is the variable used to draw, can used as center string
unit     = (/"mm/day"/)    ;this can be used as right string 
var_name = (/"series","ave"/)
case     = (/" CTRL"," NG6h"/);this can be used as left string
pre_case = (/"(a","(b"/)
years    = ispan(1979,2005,1)   ;the calculate time

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/concept_map/mdata/" + (/"CTRL","NUDG6h"/) + "-preci_4TP_annual_series.nc"
fig_name = path + "/concept_map/fig/1p5x2_preci_5area_JJAts"
fig_out  ="pdf";"x11"

area = (/"STP","CTP","NTP","SEA","NEA"/)  ;S mean south,N mean North,EA mean East Asia
lats = (/25   ,30   ,35   , 22  , 30  /) ;-20
latn = (/30   ,35   ,40   , 30  , 45  /)  ;60
lonl = (/70   ,70   ,70   ,107  ,107  /)  ;20
lonr = (/105  ,105  ,105  ,130  ,130  /);130 ;220

;area = (/"TP","STP","CTP","NTP","SEA","NEA"/)  ;S mean south,N mean North,EA mean East Asia
;lats = (/25  ,25   ,30   ,35   , 22  , 30  /) ;-20
;latn = (/40  ,30   ,35   ,40   , 30  , 45  /)  ;60
;lonl = (/70  ,70   ,70   ,70   ,107  ,107  /)  ;20
;lonr = (/105 ,105  ,105  ,105  ,130  ,130  /);130 ;220

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.01
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 5.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
ncase = dimsizes(case)
narea = dimsizes(area)
nday  = 92
nday_s=151
nday_e=242
var   = new((/ncase,narea,nyear,nday/),float)
ave   = new((/ncase,narea,      nday/),float)
plot  = new((/ncase*narea/),graphic)
plot2 = new((/ncase*narea/),graphic)
plot3 = new((/ncase*narea/),graphic)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nc = 0,ncase-1,1
    f    = addfile(filename(nc),"r")
    var(nc,:,:,:)  = f->$var_name(0)$(1:5,:,nday_s:nday_e)
    ave(nc,:,:)    = f->$var_name(1)$(1:5,nday_s:nday_e )
end do   ;the loop of the three cases

printVarSummary(var)
printVarSummary(ave)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    resxy     = True
    resxy@vpWidthF  = 0.45
    resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = ispan(4,56,2)
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    ;resxy@tmXBLabelFontHeightF     = tick_font-0.01
    ;resxy@tmXBLabelFontThicknessF  = tick_thick-1.0 
    ;resxy@tmYRLabelFontHeightF     = tick_font
    ;resxy@tmYRLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 90
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(0,95,5)
    resxy@tmXBLabels          = ispan(0,95,5)
    ;resxy@tmXBMode            = "Automatic"
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
    resxy@gsnRightStringFontHeightF  = subtitle_font
    resxy@gsnLeftStringFontHeightF   = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""

    resxy2 = resxy
    resxy2@xyMonoLineColor  = True
    resxy2@xyMonoLineThickness = True
    resxy2@xyLineThicknessF = 4.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy2@xyLineColor      = "black"

    resxy3 = resxy2
    resxy3@xyLineThicknessF = 6.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

do na =0, narea-1,1
do nc = 0, ncase-1,1
    nn = na +1
    np = ncase*na + nc
    resxy@gsnLeftString   = pre_case(nc) + nn + case(nc)
    resxy@gsnRightString  = area(na) + " (" + unit + ")"  
    resxy@trYMinF = floor(min(var(:,na,:,:)))
    resxy@trYMaxF =  ceil(max(var(:,na,:,:)))
    plot(np)  = gsn_csm_xy(wks,ispan(0,91,1),(/var(nc,na,:,:)/),resxy)
    plot2(np) = gsn_csm_xy(wks,ispan(0,91,1),(/ave(0 ,na,  :)/),resxy2)
    plot3(np) = gsn_csm_xy(wks,ispan(0,91,1),(/ave(nc,na,  :)/),resxy3)
    overlay(plot(np),plot2(np))
    overlay(plot(np),plot3(np))
end do
end do
    
   gsn_panel(wks,plot,(/narea,ncase/),resp)
end

