;***************************************************************************
;filename:      draw_1p3x2_preci_JJAts.ncl
;plot the JJA times of every year's annual cycle(id=151~242) of TP average UVT(27 curves) 
;then overlaying the climate annual cycle(black bold curve)
;and one panel has six plots.
;                                            by Ql_Ren
;                                           2017/03/29
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
draw_var = (/"preci" ,"static"  /)  ;this is the variable used to draw, can used as center string
unit     = (/"mm/day","K/100hPa"/)    ;this can be used as right string 
var_name = (/"series","ave"/)
case     = (/" CTRL"," NG6h","NG24h"/);" NG24h",this can be used as left string
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
years    = ispan(1979,2005,1)   ;the calculate time
nv = 1

;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/concept_map/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-" + draw_var(nv) + "_2EA_annual_series.nc"
fig_name = path + "/concept_map/fig/1p3x2_" + draw_var(nv) + "_2EA_JJAts"
;fig_out  = "eps";"pdf";"x11"
fig_out  = "pdf";"x11"

;area = (/"SETP","SCTP","SWTP","CETP","CWTP","NTP"/);"ESTP","ECTP","NTP","SEA","NEA"/) "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
;lats = (/25    ,25    ,25    ,30    ,30    ,35   /);25    ,30    ,35   , 22  , 30  /) 25  ,25   ,30   , ;-20
;latn = (/30    ,30    ,30    ,35    ,35    ,40   /);30    ,35    ,40   , 30  , 45  /) 40  ,30   ,35   , ;60
;lonl = (/100   ,85    ,70    ,90    ,70    ,70   /);100   ,90    ,70   ,107  ,107  /) 70  ,70   ,70   , ;20
;lonr = (/105   ,100   ,85    ,105   ,90    ,105  /);105   ,105   ,105  ,130  ,130  /) 105 ,105  ,105  , ;130 ;220
area = (/"SEA","NEA"/); (/"SETP","CETP","NTP"/);"ESTP","ECTP","NTP","SEA","NEA"/) "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 22  , 30  /); (/25    ,30    ,35   /);25    ,30    ,35   , 22  , 30  /) 25  ,25   ,30   , ;-20
latn = (/ 30  , 45  /); (/30    ,35    ,40   /);30    ,35    ,40   , 30  , 45  /) 40  ,30   ,35   , ;60
lonl = (/107  ,115  /); (/100   ,90    ,70   /);100   ,90    ,70   ,107  ,107  /) 70  ,70   ,70   , ;20
lonr = (/130  ,130  /); (/105   ,105   ,105  /);105   ,105   ,105  ,130  ,130  /) 105 ,105  ,105  , ;130 ;220

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
ncase = dimsizes(case)
narea = dimsizes(area)
nday  =121; 90
nday_s=121;152
nday_e=241
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
    var(nc,:,:,:)  = f->$var_name(0)$(:,:,nday_s:nday_e)
    ave(nc,:,:)    = f->$var_name(1)$(:,nday_s:nday_e )
end do   ;the loop of the three cases

printVarSummary(var)
printVarSummary(ave)
printVarSummary(var&day)

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
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 90
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(nday_s+1,nday_e+1,30)
    resxy@tmXBLabels          = (/"1May","31May","30Jun","30Jul","29Aug"/)
    resxy@trXMaxF      = nday_e+1;242
    resxy@trXMinF      = nday_s+1;153
    ;resxy@tmXBMode            = "Automatic"
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    if(nv.eq.1) then 
    resxy@gsnYRefLine           = 0.0
    resxy@gsnYRefLineColor      = "black"
    resxy@gsnYRefLineThicknessF = 3.0
    resxy@gsnYRefLineDashPattern= 1
    end if

    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font
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

do nc = 0, ncase-1,1
do na =0, narea-1,1
    np = nc*narea + na
    resxy@gsnLeftString   = pre_case(nc) + ")" + case(nc)
    resxy@gsnRightString  = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"  
    resxy@trYMinF = -6;floor(min(var(:,na,:,:)))
    resxy@trYMaxF =  4;ceil(max(var(:,na,:,:)))
    plot(np)  = gsn_csm_xy(wks,var&day,(/var(nc,na,:,:)/),resxy)
    plot2(np) = gsn_csm_xy(wks,var&day,(/ave(0 ,na,  :)/),resxy2)
    plot3(np) = gsn_csm_xy(wks,var&day,(/ave(nc,na,  :)/),resxy3)
    overlay(plot(np),plot2(np))
    overlay(plot(np),plot3(np))
end do
end do
    
   ;gsn_panel(wks,plot,(/narea,ncase/),resp)
   gsn_panel(wks,plot,(/ncase,narea/),resp)
end

