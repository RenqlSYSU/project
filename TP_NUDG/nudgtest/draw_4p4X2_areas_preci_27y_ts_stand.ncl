;***************************************************************************
;filename:      draw_4p4X2_areas_preci_27y_ts_stand.ncl
;input_data:    CTRL-Areas_preci_27y_ts.nc,NUDG6h-Areas_preci_27y_ts.nc 
;read the seasonly preci time series of seven areas
;and then draw the xy plot,will have 4 panels and one panel
;has seven plots with every one has three cases.
;                                            by Ql_Ren
;                                           2018/01/15
;*****************************************************************************
begin
;the path of the file
path   = "."
filein = path + (/"/mdata/CTRL-Areas_preci_27y_ts.nc","/mdata/NUDG24h-Areas_preci_27y_ts.nc","/mdata/NUDG6h-Areas_preci_27y_ts.nc"/)
fig_name = path + "/fig/4p4X2_areas_preci_27y_ts_stand"
fig_out  = "pdf"
fig_title= "Preci Time Series(stand) in " ;+ season
case     = (/"AMIP_CTRL","TP_NUDG-24h","TP_NUDG-6h"/) ; the legend name
unit     = (/"mm/day"/) ;the right title
var_name = (/"season_series"/)   ;the left title 
season   = (/"DJF","MAM","JJA","SON"/)
area  = (/"India","Tibet","Philippines","MaritimeContinent",\
          "NorthChina","CentralEastChina","SouthChina","TPnudg"/)
years    = ispan(1979,2005,1)   ;the calculate time

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 3.0
   mp_thick       = 5.0

ncase = dimsizes(case)
nseason = dimsizes(season)
narea = dimsizes(area)
nyear = dimsizes(years)
var   = new((/ncase,narea,nyear,nseason/),float)
plot  = new((/narea/),graphic)

;==============================================
;read data var(nvar,ncase,nyear)
;===================================================
do nc = 0, ncase-1,1
    f    = addfile(filein(nc),"r")
    var(nc,:,:,:)  = f->$var_name$  ;var(ncase,narea,nyear,nseason)
end do
;ave  = dim_avg_n_Wrap(var,2) ;get (ncase,narea,nseason)
;vari = dim_variance_n_Wrap(var,2) ;get (ncase,narea,nseason)
var  = dim_standardize_n_Wrap(var,1,2) ;get standardized data
printVarSummary(var)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.4
    resxy@vpHeightF = 0.2
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "landscape"
    
    resxy@xyDashPatterns     = (/0,0,0/)
    resxy@xyLineThicknesses = (/4.0, 4.0,4.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"blue","red","green"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    resxy@trYMinF = -3.0
    resxy@trYMaxF = 3.0
    
    resxy@tmXBLabelFontHeightF     = tick_font ;-0.01
    resxy@tmXBLabelFontThicknessF  = tick_thick ;-1.0 
    resxy@tmYRLabelFontHeightF     = tick_font
    resxy@tmYRLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = True ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    resxy@tmXBMode           = "Manual"
    resxy@tmXBTickSpacingF   = 5

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 4.0
    resxy@tmXMinorGrid                   = True
    resxy@tmXMinorGridLineDashPattern    = 2
    resxy@tmXMinorGridLineColor          = "gray"
    resxy@tmXMinorGridThicknessF         = 2.0
                       
    resxy@gsnYRefLine             = 0.0
    resxy@gsnYRefLineDashPatterns = 12   
    resxy@gsnYRefLineThicknesses  = 6.0
    resxy@gsnYRefLineColors       = "black"

    ;resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
    resxy@gsnRightStringFontHeightF  = subtitle_font
    resxy@gsnLeftStringFontHeightF   = subtitle_font
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.06
    resxy@pmLegendHeightF  = 0.06
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.4   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.3    ; move ref vector along X-axis

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

    resxy@gsnRightString   = ""
do ns = 0, nseason-1,1
    do na = 0, narea-1,1
        ;resxy@gsnCenterString  = "CTRL:"+vari(nv,0,ns,nl)+"  NUDG6h:"+vari(nv,1,ns,nl)+"  NUDG24h:"+vari(nv,2,ns,nl)
        resxy@gsnLeftString    = area(na)
        plot(na) = gsn_csm_xy(wks,var&year,(/var(:,na,:,ns)/),resxy)
    end do
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title + season(ns) ;+ " " + lev(nl) + "hPa"
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    resp@gsnPanelRowSpec  = False
   gsn_panel(wks,plot,(/4,2/),resp)
end do
end

