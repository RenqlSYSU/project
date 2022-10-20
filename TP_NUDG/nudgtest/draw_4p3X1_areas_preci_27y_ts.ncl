;***************************************************************************
;filename:      draw_4p3X1_areas_preci_27y_ts.ncl
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
;path   = "/users/yangsong3/renql/project/TP_NUDG/nudgtest"
path   = "."
filein = (/"/mdata/CTRL-Areas_preci_27y_ts.nc","/mdata/NUDG24h-Areas_preci_27y_ts.nc","/mdata/NUDG6h-Areas_preci_27y_ts.nc"/)
fig_name = path + "/fig/4p3X1_areas_preci_27y_ts"
fig_out  = "pdf"
fig_title= "Preci Time Series in " ; + season

output_write = 0 ;1 means write some data to fname txt file
fname = path + "/areas_seasonly_preci_ave_vari.txt"

case     = (/"AMIP_CTRL","TP_NUDG-24h","TP_NUDG-6h"/) ; the legend name
unit     = (/"mm/day"/) ;the right title
var_name = (/"season_series"/)   ;the left title 
season   = (/"DJF","MAM","JJA","SON"/)
area  = (/"India","Tibet","Philippines","MaritimeContinent",\
          "NorthChina","CentralEastChina","SouthChina"/)
years    = ispan(1979,2005,1)   ;the calculate time

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

ncase = dimsizes(case)
nseason = dimsizes(season)
narea = dimsizes(area)
nyear = dimsizes(years)
var   = new((/ncase,narea,nyear,nseason/),float)
plot  = new((/ncase/),graphic)

;==============================================
;read data var(narea,nyear,nseason)
;===================================================
do nc = 0, ncase-1,1
    f    = addfile(path + filein(nc),"r")
    var(nc,:,:,:)  = f->$var_name$  ;var(ncase,narea,nyear,nseason)
end do
ave  = dim_avg_n_Wrap(var,2) ;get (ncase,narea,nseason)
vari = dim_variance_n_Wrap(var,2) ;get (ncase,narea,nseason)
printVarSummary(var)
printVarSummary(vari)

;----write table to txt file----------
if(output_write .eq. 1) then
write_table(fname,"w",[/"DIF(ave,vari),MAM(ave,vari),JJA(ave,vari),SON(ave,vari)"/],"%s")
do na = 0, narea-1,1
    write_table(fname,"a",[/"  "/],"%s")
    alist = [/area(na),case,ave(:,na,0),vari(:,na,0),\
                            ave(:,na,1),vari(:,na,1),\
                            ave(:,na,2),vari(:,na,2),\
                            ave(:,na,3),vari(:,na,3)/]
    write_table(fname,"a",alist,"%s %s %10.6f%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f")
end do
end if
;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.3
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
do nc = 0, ncase-1,1
    resxy@xyExplicitLegendLabels = area + " " + ave(nc,:,ns) + " "+ vari(nc,:,ns)
    resxy@gsnLeftString    = case(nc)
    plot(nc) = gsn_csm_xy(wks,years,(/var(nc,:,:,ns)/),resxy)
end do 
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title + season(ns) 
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
   gsn_panel(wks,plot,(/3,1/),resp)
end do
end

