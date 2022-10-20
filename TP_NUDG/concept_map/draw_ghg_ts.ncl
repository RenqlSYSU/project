;***************************************************************************
;filename:      draw_ghg_ts.ncl_
;draw the curves of the development of five greenhouse gases from 1979 
;to 2005 in one plot
;                                            by Ql_Ren
;                                           2017/04/08
;******************************************************************************
begin
;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = "/users/yangsong3/CESM/input/atm/cam/ggas/ghg_hist_1765-2005_c091218.nc" 
fig_name   = path + "/concept_map/fig/ghg_ts"
fig_out    ="pdf";"x11"

years    = ispan(1979,2005,1)   ;the calculate time
draw_var = (/"CO~B~2~N~","CH~B~4~N~","N~B~2~N~O","CFC-11","CFC-12"/)
unit     = (/"10~S~-6~N~","10~S~-9~N~","10~S~-9~N~","10~S~-12~N~","10~S~-12~N~"/)+" mol/mol"  

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
   font_index     = 22

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
nyear = dimsizes(years)
nvar  = dimsizes(draw_var)
var   = new((/nvar,nyear/),float)
plot  = new((/1/),graphic)

f        = addfile(filename,"r")
var_name = getfilevarnames(f)
print(var_name)

date = f->date
date = date/10000    ;convert yyyymmdd to yyyy for date's type is integer
select_time = ind(date.ge.1979.and.date.le.2005)
printVarSummary(select_time)

do nv = 0,nvar-1,1
nvv = nv + 2
var(nv,:)  = f->$var_name(nvv)$(select_time)     ;(nv,:,:,nl)
end do
printVarSummary(var)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)
    
    gsn_define_colormap(wks,"matlab_jet") 
     ;has 64 colors blue green yellow,red

    resxy     = True
    ;resxy@vpWidthF  = 0.9
    ;resxy@vpHeightF = 0.25
    ;resxy@gsnFrame      = False
    ;resxy@gsnDraw       = False
    ;resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "landscape";"portrait"
    
    resxy@xyDashPatterns    = (/0  ,0  ,0  ,0  ,0  /)
    resxy@xyLineThicknesses = (/5.0,5.0,5.0,5.0,5.0/)   ; make second line thicker
    resxy@xyLineColors      = (/2  , 17, 33, 49, 57/)
    resxy@tiXAxisOn = False
    resxy@tiXAxisFont   = font_index
    resxy@tiXAxisString = "Year"
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))

    resxy@tmBorderThicknessF        = mp_thick  
    resxy@tmXBLabelFont             = font_index
    resxy@tmYLLabelFont             = font_index
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
    resxy@tmXBTickStartF     = 1979
    resxy@tmXBTickEndF       = 2005
    ;resxy@tmXBMode            = "Explicit"
    ;resxy@tmXBValues          = ispan(1979,2005,3)
    ;resxy@tmXBLabels          = ispan(1979,2005,3)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = draw_var + " " + unit
    resxy@lgLabelFont            = font_index 
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.15
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = -0.05   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = -0.05    ; move ref vector along X-axis

    resxy@gsnStringFont              = font_index
    resxy@gsnStringFontHeightF       = subtitle_font 
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""
    ;resxy@trYMinF = floor(min(var(:,:,:,nl)))
    ;resxy@trYMaxF =  ceil(max(var(:,:,:,nl)))
    plot  = gsn_csm_xy(wks,date(select_time),var,resxy)
    
    ;resp  = True    
    ;resp@gsnMaximize   = True
    ;resp@txString      = fig_title + lev(nl) + "hPa " + draw_var(nv)
    ;resp@txFontHeightF = title_font 
    ;resp@txFontThicknessF = title_thick
   ;gsn_panel(wks,plot,(/3,1/),resp)
end

