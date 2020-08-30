;***************************************************************************
;filename:      draw_12p3X1_27y_ts.ncl
;input_data:    Areas_preci_27y_ts.nc 
;read the seasonly U,V,T time series at three levs
;and then draw the xy plot,will have 12 panels. 
;                                            by Ql_Ren
;                                           2018/01/08
;*****************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 5.0

;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
case   = "ERA-Interim"
var_name = (/"uppEKE-2.5-10day","uppEKE-10-90day","uwnd","hgt"/)   ;,"Q1"   the left title 
unit     = (/"10~S~4~N~J"      ,"10~S~4~N~J"     ,"m/s" ,"gpm"/)   ;,"W/kg"
draw_var = (/"STEhgh"          ,"STEhgh"         ,"Jet" ,"NBR"/)   ;North Ridge; TPH
stand = True
nv = 0

;area = all, TP, upstream
lats = (/25 ,25 ,25 ,/);-20
latn = (/40 ,40 ,40 ,/);60
lonl = (/35 ,70 ,35 ,/);20
lonr = (/105,105,70 ,/);220

filename = path + "mdata/" + case + "_Clim_month_ave_" + var_name + ".nc"
if(stand) then 
fig_name = path + "fig/STE_indices_annual_ts_stand"
else
fig_name = path + "fig/STE_indices_annual_ts"
end if
fig_out  = "pdf"

nvar = dimsizes(var_name)
ts   = new((/nvar,12/),float)
plot = new((/nvar/),graphic)
;==============================================
;read data var(nvar,ncase,nyear)
;===================================================
nv = 0
f   = addfile(filename(nv),"r")
var = f->month_ave(:,
ts  = ts/100000

if(stand) then 
var = dim_standardize_n(var,1,1)
end if
ave  = dim_avg_n_Wrap(var,1) ;ave(nvar,ncase,nseason,nlev)
vari = dim_variance_n_Wrap(var,0) ;vari(nvar,ncase,nseason,nlev)
printVarSummary(var)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.25
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns     = (/0,0,0/)
    resxy@xyLineThicknesses = (/2.0,2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"blue","red","black"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
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
    resxy@tmXBTickSpacingF   = 5

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    resxy@gsnYRefLine            = ave
    resxy@gsnYRefLineDashPatterns = (/12,12,12/)   
    resxy@gsnYRefLineThicknesses  = (/3.0,3.0,3.0/)
    resxy@gsnYRefLineColors       = (/"blue","red","black"/)
                             
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = (/"DJF","JJA","All"/) 
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.1
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.35   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.4    ; move ref vector along X-axis

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
;   
;    resxy@gsnCenterString  = vari(nv,0,ns,nl)+"  NUDG6h:"+vari(nv,1,ns,nl)+"  NUDG24h:"+vari(nv,2,ns,nl)
    resxy@gsnLeftString    = case 
    resxy@gsnRightString   = var_name(nv) + " (" + unit(nv) + ")"
    if(stand) then 
    resxy@gsnYRefLine      = (/-1.5,0.0,1.5/)
    else
    resxy@gsnYRefLine      = (/ave(0),ave(1),ave(2)/)
    end if
    plot = gsn_csm_xy(wks,years,(/var/),resxy)

    gsn_panel(wks,plot,(/1,1/),resp)
end

