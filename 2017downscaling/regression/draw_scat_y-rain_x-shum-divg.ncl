;*********************************************************************
;filename:      draw_scat_y-rain_x-var.ncl_
;read rain of Ruyuan and area average omega from 1970-2016 in winter
;then calculate the regress the rain by the omega;
;then draw the scatter and regline plot
;                                            by Ql_Ren
;                                           2018/03/03
;*********************************************************************************
begin
path     = "/home/ys17-19/renql/project/2017downscaling"
var_name = (/"rain","uwnd","vwnd","shum","pres.sfc"/)
filein1  = path + "/data/rain7016.nc"
filein2  = path + "/data/NECP_NCAR-daily/" + var_name(1:4) + "/" + var_name(1:4) + "*"
years    = ispan(1970,2016,1)   ;the calculate time
season   = "Winter"
draw_var = (/"shum-divg","Moisture Flux Divergence"/)
unit     = "g/(m*m*s)"
g        = 9.8 ;m/(s*s)

;area of field 
lats = 20  ;22  ;-15
latn = 25  ;27  ;55
lonl = 108 ;110 ;30
lonr = 113 ;118 ;150

fig_name   = path + "/regression/fig/" + "1p_scat_y-rain_x-" + draw_var(0) + "_" + season ;-delay" 
fig_title  = "Regression with " + draw_var(1) + " in " + season
fig_out    = "pdf";"x11"

;==========================================================================
;read the winter(9,10,11,12,1,2) omega data,2016-1948=68,1970-1948=22 and ruyuan rainfall
;===========================================================================
files       = systemfunc("ls " + filein2(0)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
time        = cd_calendar(f2[:]->time,0)
select_time = ind(time(:,1).ge.9.or.time(:,1).le.2)
ntime       = dimsizes(select_time)
var_u      := f2[:]->$var_name(1)$(select_time,0:7,{lats:latn},{lonl:lonr})

files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v       := f2[:]->$var_name(2)$(select_time,0:7,{lats:latn},{lonl:lonr})

files       = systemfunc("ls " + filein2(2)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var        := f2[:]->$var_name(3)$(select_time,:,{lats:latn},{lonl:lonr})
var_shum    = var*1000      ;convert unit from kg/kg to g/kg
copy_VarMeta(var,var_shum)

files       = systemfunc("ls " + filein2(3)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var        := f2[:]->pres(select_time,{lats:latn},{lonl:lonr})
pres_sfc    = var/100       ;convert unit from Pa to hPa
copy_VarMeta(var,pres_sfc)

;var2    = wgt_vert_avg_beta ( var_shum&level, var_shum, pres_sfc, 0, 0 ) / g 
;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum

var_u    = var_u*var_shum/g
var_v    = var_v*var_shum/g
var_divg = uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 3) ;calculate the divergence
var2     = wgt_vert_avg_beta ( var_shum&level, var_divg, pres_sfc, 0, 0 )  

var2 := wgt_areaave_Wrap(var2,1,1,0)

f1      = addfile(filein1,"r")
var1    = f1->$var_name(0)$(select_time)  ;get var1(time)

rc    = regCoef_n(var2(0:ntime-2),var1(1:ntime-1),0,0)
var1r = rc*var2(0:ntime-2)+rc@yintercept
prob  = student_t(rc@tval,rc@nptxy-2)
print(rc)

;====================================================================
;draw the scatter plot and regline
;=========================================================================
lefttitle   = "t-value:" + rc@tval 
righttitle  = lats + "-" + latn + "N, " + lonl + "-" + lonr +"E"
centertitle = "prob:" + prob
xaxis      = draw_var(1) + " " + unit
yaxis      = "Ruyuan daily rainfall (mm)"
yref       = 50 ;mm
legendlabel= "y = " + rc + "*x + " + rc@yintercept

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

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    ;res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "landscape"

;-------set the attribute of X Y axis----------
    res@tiXAxisOn           = True
    res@tiXAxisFontHeightF  = label_font
    res@tiXAxisFontThicknessF = 1
    res@tiXAxisFontColor    = "black"
    res@tiXAxisString       = xaxis
    
    res@tiYAxisOn           = True
    res@tiYAxisString       = yaxis

    res@trYMinF = 0
    res@trYMaxF = 120
    
    ;res@tmXBLabelFontHeightF     = tick_font-0.01
    ;res@tmXBLabelFontThicknessF  = tick_thick-1.0 
    ;res@tmYRLabelFontHeightF     = tick_font
    ;res@tmYRLabelFontThicknessF  = tick_thick

    res@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    res@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    res@tmXTOn               = False ; Turn off x-axis (top) major ticks
    res@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    res@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    res@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;res@tmXBMajorLengthF = 0.01
    ;res@tmYRMajorLengthF = 0.01
    ;res@tmXBMode           = "Manual"
    ;res@tmXBTickSpacingF   = 3

    res@tmXMajorGrid                   = True
    res@tmXMajorGridLineDashPattern    = 2
    res@tmXMajorGridLineColor          = "gray"
    res@tmXMajorGridThicknessF         = 2.0
                       
;-------set the attribute of the mark and line----------
    res@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
    res@xyMarkers           = 16                     ; choose type of marker 
    res@xyMarkerColor       = "red"                  ; Marker color
    res@xyMarkerSizeF       = 0.005                  ; Marker size (default 0.01)
    res@xyDashPatterns      = 1                      ; solid line 
    res@xyLineThicknesses   = (/1,2/)                ; set second line to 2
    res@xyLineColors        = "blue"                 ; change line color
    
;------------legend---------------------
    res@pmLegendDisplayMode = "Always"
    res@lgOrientation       = "Vertical"
    res@lgPerimOn           = False
    res@lgPerimFill         = "SolidFill"
    res@lgPerimFillColor    = "white"
    res@xyExplicitLegendLabels = (/"",legendlabel/)
    res@lgLabelFontColor       = "black"
    res@lgLabelFontThicknessF  = label_thick
    res@lgLabelFontHeightF     = label_font
    res@pmLegendWidthF   = 0.07
    res@pmLegendHeightF  = 0.1
    res@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    res@pmLegendOrthogonalPosF  = -0.35   ; move ref vector along Y-axis
    res@pmLegendParallelPosF    = 0.2    ; move ref vector along X-axis
    
    res@gsnYRefLine             = yref
    res@gsnYRefLineDashPatterns = 12   
    res@gsnYRefLineThicknesses  = 3.0
    res@gsnYRefLineColors       = "green"
                             
    res@gsnCenterStringFontHeightF = subtitle_font-0.005
    res@gsnRightStringFontHeightF  = subtitle_font-0.005
    res@gsnLeftStringFontHeightF   = subtitle_font-0.005
    
    res@gsnCenterString  = centertitle
    res@gsnLeftString    = lefttitle 
    res@gsnRightString   = righttitle
  
    plot  = gsn_csm_xy (wks,var2(0:ntime-2),(/var1(1:ntime-1),var1r/),res)        ; create plot

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    resp@txString      = fig_title
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
   gsn_panel(wks,plot,(/1,1/),resp)
end
