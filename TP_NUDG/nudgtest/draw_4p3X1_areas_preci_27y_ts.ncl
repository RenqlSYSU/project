;***************************************************************************
;filename:      draw_4p2X1_areas_preci_27y_ts.ncl
;input_data:    CTRL-Areas_preci_27y_ts.nc 
;               NUDG6h-Areas_preci_27y_ts.nc
;read the seasonly preci time series of seven areas
;and then draw the xy plot,will have 4 panels and one panel
;has two plot which 1 is CTRL lines and 2 is NUDG lines. 
;                                            by Ql_Ren
;                                           2018/01/11
;*****************************************************************************
begin
;the path of the file
path   = "/users/yangsong3/renql/project/TP_NUDG/nudgtest"
filein = (/"/mdata/CTRL-Areas_preci_27y_ts.nc","/mdata/NUDG6h-Areas_preci_27y_ts.nc"/)
fig_name = path + "/fig/4p2X1_areas_preci_27y_ts"
fig_out  = "pdf"
fig_title= "Preci Time Series in " ; + season
case     = (/"AMIP_CTRL","TP_NUDG-6h"/) ; the legend name
unit     = (/"mm/day"/) ;the right title
var_name = (/"season_series"/)   ;the left title 
season   = (/"DJF","MAM","JJA","SON"/)
area  = (/"India","Tibet","Philippines","Maritime_continent",\
          "North China","Central East China","South China"/)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 5.0

plot  = new((/2/),graphic)
ncase = dimsizes(case)
nseason = dimsizes(season)
narea = dimsizes(area)

;==============================================
;read data var(narea,nyear,nseason)
;===================================================
f    = addfile(path + filein(0),"r")
var1 = f->$var_name$  ;var(narea,nyear,nseason)
f    = addfile(path + filein(1),"r")
var2 = f->$var_name$  ;var(narea,nyear,nseason)
ave1= dim_avg_n_Wrap(var1,1) ;get (narea,nseason)
ave2= dim_avg_n_Wrap(var2,1)
vari1= dim_variance_n_Wrap(var1,1) ;get (narea,nseason)
vari2= dim_variance_n_Wrap(var2,1)
printVarSummary(var1)
printVarSummary(vari2)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.4
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns     = (/0,0,0,0,16,2,0/)
    resxy@xyLineThicknesses = (/3.0, 3.0, 3.0, 3.0, 6.0, 6.0, 6.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"blue","red","forestgreen","darkgoldenrod","black","black","black"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    resxy@trXMinF = 1979
    resxy@trXMaxF = 2005
    ;resxy@trYMinF = -3.0
    ;resxy@trYMaxF = 15.0
    
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
    resxy@tmXBMode           = "Manual"
    resxy@tmXBTickSpacingF   = 3

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    ;resxy@gsnYRefLine            = ave
    ;resxy@gsnYRefLineDashPatterns = (/12,12,12/)   
    ;resxy@gsnYRefLineThicknesses  = (/3.0,3.0,3.0/)
    ;resxy@gsnYRefLineColors       = (/"blue","red","green"/)
                             
    resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
    resxy@gsnRightStringFontHeightF  = subtitle_font
    resxy@gsnLeftStringFontHeightF   = subtitle_font
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    ;resxy@xyExplicitLegendLabels = area
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.09
    resxy@pmLegendHeightF  = 0.3
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.8    ; move ref vector along X-axis

; Create legend
;    lgres                    = True
;    lgres@lgLineColors       = (/"blue","red"/)
;    lgres@lgLineThicknessF   = (/2.0,   2.0/) 
;    lgres@lgLabelFontHeightF = 0.1            ; set the legend label font thickness
;    lgres@vpWidthF           = 0.1           ; width of legend (NDC)
;    lgres@vpHeightF          = 0.05            ; height of legend (NDC)
;    lgres@lgMonoDashIndex    = True   
;    lgres@lgPerimOn          = True          ; draw the box perimeter in orange
;    lgres@lgLabelJust = "TopLeft"
;    labels = " " + case

    resxy@gsnRightString   = unit
do ns = 0, nseason-1,1
    resxy@xyExplicitLegendLabels = area + " " + ave1(:,ns) + " "+ vari1(:,ns)
    resxy@gsnLeftString    = case(0)
    plot(0) = gsn_csm_xy(wks,var1&year,(/var1(:,:,ns)/),resxy)
    resxy@xyExplicitLegendLabels = area + " " + ave2(:,ns) + " "+ vari2(:,ns)
    resxy@gsnLeftString    = case(1)
    plot(1) = gsn_csm_xy(wks,var2&year,(/var2(:,:,ns)/),resxy)
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title + season(ns) 
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
   gsn_panel(wks,plot,(/2,1/),resp)
end do
end

