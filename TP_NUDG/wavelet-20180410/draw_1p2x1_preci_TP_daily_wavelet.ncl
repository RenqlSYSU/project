;*********************************************************************************
;filename:      draw_1p2x1_preci_TP_daily_wavelet.ncl
;input_data:    CTRL-preci_TP_JJAdaily_wavelet.nc
;                                            by Ql_Ren
;                                           2018/04/10
;*********************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.015
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

;path       = "/home/ys17-19/renql/project/TP_NUDG"
path       = "/users/yangsong3/renql/project/TP_NUDG"
;filename   = path + "/wavelet-20180410/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-preci_TP_JJAdaily_wavelet.nc"
;fig_name   = path + "/wavelet-20180410/fig/1p3x1_preci_TP_daily_wavelet_nobias"
      ;case = (/"(a) CTRL","(b) NG24h","(c) NG6h"/)
filename   = path + "/wavelet-20180410/mdata/" + (/"CTRL","NUDG6h"/) + "-preci_TP_JJAdaily_wavelet.nc"
fig_name   = path + "/wavelet-20180410/fig/1p2x1_preci_TP_daily_wavelet_nobias"
      case = (/"(a) CTRL","(b) NG6h"/)
fig_out    = "pdf" ;"eps";"x11"
draw_var   = "preci"
unit       = "mm~S~2~N~/day~S~2~N~"
var_name   = (/"wave"/)

;area:   TP , east Asia
lats = (/25 , 20/) ;-20
latn = (/40 , 45/)  ;60
lonl = (/70 ,107/)  ;20
lonr = (/105,122/);130 ;220
na = 0

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
npart  = getfilevardimsizes(f,"part")
njtot  = getfilevardimsizes(f,"jtot")
ntime  = getfilevardimsizes(f,"time")
power  = new((/ncase,njtot,ntime/),float)
sig    = new((/ncase,njtot,ntime/),float)
plot   = new((/ncase/),graphic)
plot2  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
do nf = 0,ncase-1,1
    f   = addfile(filename(nf),"r")
    var = f->$var_name$       ;(ncase,neval,nlat,nlon)
    power(nf,:,:) = onedtond(var@power,(/njtot,ntime/))   
      sig(nf,:,:) = power(nf,:,:)/conform(power(nf,:,:),var@signif,0) ;compute significance ( >= 1 is significant)
    power(nf,:,:) = power(nf,:,:)/conform(power(nf,:,:),var@scale ,0)   ;no bias, referenced by Liu et.al(2007)
end do
power!1         = "period"
power&period    = var@period
power!2         = "time"
power&time      = var&time
copy_VarMeta(power,sig)
sig = mask(sig,sig.ge.1.0,True)
printVarSummary(power)
printVarSummary(sig)

exponent = ispan(2,11,2)
ylabel   = exp( exponent*log(2.0) )
print(ylabel)
;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@vpHeightF        = 0.2                   ; height and width of plot
    res@vpWidthF         = 0.9

;------------set the form of the coordinate-----------------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXTMajorOutwardLengthF = 0.01
    res@tmXBMajorOutwardLengthF = 0.01
    res@tmYLMajorOutwardLengthF = 0.01
    res@tmYRMajorOutwardLengthF = 0.01
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    res@tmXTMinorOutwardLengthF = 0.01
    res@tmXBMinorOutwardLengthF = 0.01
    res@tmYLMinorOutwardLengthF = 0.01
    res@tmYRMinorOutwardLengthF = 0.01

;---------set the attribute of the x y axis------------
    res@tiXAxisOn       = False
    res@tiXAxisString   = "Year"
    
    res@tiYAxisOn       = True
    res@tiYAxisFont     = font_index
    res@tiYAxisFontHeightF = subtitle_font 
    res@tiYAxisString   = "period(days)"
    
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = ispan(1,ntime,92)
    res@tmXBLabels      = ispan(1979,2005,1)
    
    res@trYReverse          = True   ; reverse y-axis
    res@tmYLMode            = "Explicit"
    res@tmYLValues          = ylabel 
    res@tmYLLabels          = ylabel
    res@tmLabelAutoStride   = True

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    ;res@cnFillPalette       = "BlAqGrYeOrReVi200"   ; set color map
    ;res@cnFillMode          = "RasterFill"          ; turn on raster mode
    ;res@cnRasterSmoothingOn = True                  ; turn on raster smoothing
    res@cnLineLabelsOn   = False  ; True
    res@cnFillMode          = "RasterFill"          ; turn on raster mode
    res@cnRasterSmoothingOn = True                  ; turn on raster smoothing

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    res@cnLevels             = (/0.01,0.03,0.06,0.1 ,0.15, \
                                 0.2 ,0.3 ,0.4 ,0.6 ,0.8 , \
                                 1   ,1.2 ,1.5 ,2   ,2.5 , 3 , 4/)
    ;res@cnLevels             = ispan(0,60,10)
    colors       = new((/18/),integer)
    colors(0:8)  = ispan(2,10,1)
    colors(9:17) = ispan(13,21,1)
    res@cnFillColors         = colors 
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = -0.09       ; set min contour level
    ;res@cnMaxLevelValF        = 0.09       ; set max contour level
    ;res@cnLevelSpacingF       = 0.01        ; set contour interval

;    res@gsnSpreadColors      = True              ; turn of the full color map
;    res@gsnSpreadColorStart  = 2                ; start at color 17
;    res@gsnSpreadColorEnd    = 21                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"   ;"Vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-2
    
;----------Set res sig-----------
    res_sig                  = True
    res_sig@gsnDraw          = False
    res_sig@gsnFrame         = False
    
    res_sig@tiXAxisOn       = False
    res_sig@tiYAxisOn       = False
    res_sig@tmXBMode        = "Explicit"
    res_sig@tmXBValues      = ispan(1,ntime,92)
    res_sig@tmXBLabels      = ispan(1979,2005,1)
    res_sig@trYReverse          = True   ; reverse y-axis
    res_sig@tmYLMode            = "Explicit"
    res_sig@tmYLValues          = ylabel
    res_sig@tmYLLabels          = ylabel
    res_sig@tmLabelAutoStride   = True

    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "grey0"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    
    res_sig@gsnCenterString      = ""
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString   = draw_var + " (" + unit +")"

    resp  = True    
    resp@gsnMaximize   = True
;    resp@txFont        = font_index
;    resp@txFontHeightF = title_font 
;    resp@txFontThicknessF = title_thick
;    resp@txString      =  fig_title 
    resp@gsnPaperOrientation = "landscape"

do nc = 0,ncase-1,1
    res@gsnLeftString    = case(nc)
    plot(nc)  = gsn_csm_contour(wks,power(nc,:,:),res)
    plot2(nc) = gsn_csm_contour(wks,sig(nc,:,:)  ,res_sig)
    overlay(plot(nc),plot2(nc))
    print(max(power(nc,:,:)))
    print(min(power(nc,:,:)))
end do

gsn_panel(wks,plot,(/ncase,1/),resp)
end

