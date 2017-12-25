;***************************************************************************
;read the seasonly U,V,T time series at three levs
;and then draw the xy plot,will have 12 panels. 
;                                            by Ql_Ren
;                                           2017/12/07
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
filein   = "/users/yangsong3/Model_Team/F/nudgtest/mdata/27years_time_series.nc"
fig_name = "/users/yangsong3/Model_Team/F/nudgtest/fig/12panel_3X1_27years_time_series"
fig_out  = "pdf"
fig_title= "Time series of TP area at " ; + season + lev
case     = (/"AMIP_CTRL","TP_NUDG"/) ; the legend name
unit     = (/"m/s","m/s","~S~o~N~C"/) ;the right title
var_name = (/"U","V","T"/)   ;the left title 
season   = (/"DJF","MAM","JJA","SON"/)
lev      = (/850,500,200/)

plot  = new((/3/),graphic)
nvar  = dimsizes(var_name)
ncase = dimsizes(case)
nseason = dimsizes(season)
nlev  = dimsizes(lev)

;==============================================
;read data var(nvar,ncase,nyear)
;===================================================
f   = addfile(filein,"r")
var = f->season_series  ;var(nvar,ncase,nyear,nseason,nlev)
var(2,:,:,:,:) = var(2,:,:,:,:)-273
ave = dim_avg_n_Wrap(var,2) ;ave(nvar,ncase,nseason,nlev)
vari= dim_variance_n_Wrap(var,2) ;vari(nvar,ncase,nseason,nlev)
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
    
    resxy@xyDashPatterns     = (/0,0/)
    resxy@xyLineThicknesses = (/2.0, 2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"blue","red"/) ; change line color
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
    resxy@tmXBTickSpacingF   = 3

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    ;resxy@gsnYRefLine            = ave
    resxy@gsnYRefLineDashPatterns = (/12,12/)   
    resxy@gsnYRefLineThicknesses  = (/3.0,3.0/)
    resxy@gsnYRefLineColors       = (/"blue","red"/)
                             
    resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
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
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.07
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.38   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.4    ; move ref vector along X-axis

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

do ns = 0, nseason-1,1
do nl = 0, nlev-1,1
    do nv = 0, nvar-1,1
        resxy@gsnCenterString  = "CTRL:"+vari(nv,0,ns,nl)+"  NUDG:"+vari(nv,1,ns,nl)
        resxy@gsnLeftString    = var_name(nv)
        resxy@gsnRightString   = unit(nv)
        resxy@gsnYRefLine      = (/ave(nv,0,ns,nl),ave(nv,1,ns,nl)/)
        plot(nv) = gsn_csm_xy(wks,var&year,(/var(nv,:,:,ns,nl)/),resxy)
    end do
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title + season(ns) + " " + lev(nl) + "hPa"
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
   gsn_panel(wks,plot,(/3,1/),resp)
end do
end do
end

