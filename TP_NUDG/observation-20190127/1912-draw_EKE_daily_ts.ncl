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
var_name = (/"EKEts"/)   ;the left title 
unit     = (/"10~S~5~N~J"/)
stand = True
years = ispan(1979,2010,1)
nv = 0

filename = path + "mdata/" + case + "_intEKE_TP_ts-2.5-10day.nc"
if(stand) then 
fig_name = path + "fig/EKE_TP_ts-stand"
else
fig_name = path + "fig/EKE_TP_ts"
end if
fig_out  = "pdf"
;fig_title= "Time series of TP area at " ; + season + lev

nday = 30
ms = 5 ;start month
nyear = dimsizes(years)
var = new((/4,nyear*nday/),float)
plot  = new((/4/),graphic)

;==============================================
;read data var(nvar,ncase,nyear)
;===================================================
f   = addfile(filename,"r")
time = cd_calendar(f->time,0)
do nm=0,3,1
    nm1 = nm + ms
    select_time := ind(time(:,1).eq.5)
    ts := f->$var_name(nv)$(0:nyear-1,select_time)
    ts = dim_standardize_n(ts,1,0)
    var(nm,:) = ndtooned(ts(:,0:(nday-1)))
end do
;var = var/100000

if(stand) then 
;var = dim_standardize_n(var,1,1)
end if
ave  = dim_avg_n_Wrap(var,1) ;ave(nvar,ncase,nseason,nlev)
vari = dim_variance_n_Wrap(var,1) ;vari(nvar,ncase,nseason,nlev)
printVarSummary(var)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.2
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns     = (/0/)
    resxy@xyLineThicknesses = (/2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black"/) ; change line color
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
    resxy@tmXBTickSpacingF   = 90

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    ;resxy@gsnYRefLine            = ave
    resxy@gsnYRefLineDashPatterns = (/12,12,12/)   
    resxy@gsnYRefLineThicknesses  = (/3.0,3.0,3.0/)
    resxy@gsnYRefLineColors       = (/"red","black","red"/)
                             
    ;legend
    resxy@pmLegendDisplayMode = "NoCreate" ;"Always"
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
    if(stand) then 
    resxy@gsnYRefLine      = (/-1.5,0.0,1.5/)
    else
    resxy@gsnYRefLine      = (/ave(0),ave(1),ave(2)/)
    end if

x = ispan(1,nyear*nday,1)
do np = 0, 3, 1
    nm1 = np + ms
    resxy@gsnRightString   = nm1 + " " + var_name(nv) + " (" + unit(nv) + ")"
    plot(np) = gsn_csm_xy(wks,x,(/var(np,:)/),resxy)
end do
    gsn_panel(wks,plot,(/4,1/),resp)
end

