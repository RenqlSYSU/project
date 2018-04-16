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

area =(/"STP","CTP","NTP","SEA","NEA"/)  ;S mean south,N mean North,EA mean East Asia
lats = (/25  ,30   ,35   , 22  , 30  /) ;-20
latn = (/30  ,35   ,40   , 30  , 45  /)  ;60
lonl = (/70  ,70   ,70   ,107  ,107  /)  ;20
lonr = (/105 ,105  ,105  ,130  ,130  /);130 ;220
na = 0;5;4;3;2;1
isnormal = 1 ;1 mean normalized
narea = dimsizes(area)

if (isnormal.eq.1) then
cnlevels =  (/0.002,0.005 ,0.009,0.015,0.02 , \
              0.03 ,0.05  ,0.08 ,0.12 ,0.18 , \
              0.25 ,0.35  , 0.5 , 0.7 ,1    , 1.3,1.6/)
else
cnlevels  = (/0.01,0.03,0.06,0.1 ,0.15, \
              0.2 ,0.3 ,0.4 ,0.6 ,0.8 , \
              1   ,1.2 ,1.5 ,2   ,2.5 , 3 , 4/)
end if

path       = "/users/yangsong3/renql/project/TP_NUDG"
filename   = path + "/wavelet-20180410/mdata/" + (/"CTRL","NUDG6h"/) + "-preci_" + area(na) + "_JJAdaily_wavelet.nc"
if(isnormal.eq.1) then
fig_name   = path + "/wavelet-20180410/fig/1p5x2_preci_daily_wavelet_nobias-stand"
else
fig_name   = path + "/wavelet-20180410/fig/1p5x2_preci_daily_wavelet_nobias"
end if
  pre_case = (/"(a","(b"/)  
      case = (/"CTRL","NG6h"/)  
fig_out    = "pdf" ;"eps";"x11"
draw_var   = "preci"
unit       = "mm~S~2~N~/day~S~2~N~"
var_name   = (/"wave"/)

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
npart  = getfilevardimsizes(f,"part")
njtot  = getfilevardimsizes(f,"jtot")
ntime  = getfilevardimsizes(f,"time")
power  = new((/narea,ncase,njtot,ntime/),float)
sig    = new((/narea,ncase,njtot,ntime/),float)
plot   = new((/narea*ncase/),graphic)
plot2  = new((/narea*ncase/),graphic)

;================================================================
;read data
;====================================================================
do na = 0,narea-1,1 
filename   = path + "/wavelet-20180410/mdata/" + (/"CTRL","NUDG6h"/) + "-preci_" + area(na) + "_JJAdaily_wavelet.nc"
do nf = 0,ncase-1,1
    f   = addfile(filename(nf),"r")
    var = f->$var_name$       ;(ncase,neval,nlat,nlon)
    power(na,nf,:,:) = onedtond(var@power,(/njtot,ntime/))   
      sig(na,nf,:,:) = power(na,nf,:,:)/conform(power(na,nf,:,:),var@signif,0) ;compute significance ( >= 1 is significant)
    power(na,nf,:,:) = power(na,nf,:,:)/conform(power(na,nf,:,:),var@scale ,0)   ;no bias, referenced by Liu et.al(2007)
    if(isnormal.eq.1) then
        power(na,nf,:,:) = power(na,nf,:,:)/var@stdev/var@stdev
    end if
end do
end do

power!2         = "period"
power&period    = var@period
power!3         = "time"
power&time      = var&time
copy_VarMeta(power,sig)
sig = mask(sig,sig.ge.1.0,True)
printVarSummary(power)
printVarSummary(sig)
printVarSummary(var)

exponent = ispan(2,9,2)
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
    res@vpHeightF        = 0.15                   ; height and width of plot
    res@vpWidthF         = 0.45

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
    
    res@tiYAxisOn       = False
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

    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "black"
    res@tmYMajorGridThicknessF         = 4.0

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    ;res@cnFillPalette       = "BlAqGrYeOrReVi200"   ; set color map
    ;res@cnFillMode          = "RasterFill"          ; turn on raster mode
    ;res@cnRasterSmoothingOn = True                  ; turn on raster smoothing
    res@cnLineLabelsOn   = False  ; True
    res@cnInfoLabelOn    = False
    res@cnFillMode          = "RasterFill"          ; turn on raster mode
    res@cnRasterSmoothingOn = True                  ; turn on raster smoothing

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    res@cnLevels             = cnlevels
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

    res@lbLabelBarOn         = False  ; the label bar
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

    res_sig@cnFillOn         = True;False       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnInfoLabelOn    = False
;    res_sig@cnLineThicknessF = mp_thick
;    res_sig@cnLineColors     = "black"
;    res_sig@cnLevelSelectionMode = "ExplicitLevels"
;    res_sig@cnLevels         = 
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "grey0"
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    
    res_sig@gsnCenterString      = ""
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;------------set the attribute of the cone of influence---------------------
    rescoi              = True
    rescoi@gsFillIndex  = 1 

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnCenterString   = "" ;draw_var + " (" + unit +")"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Horizontal" ;"Vertical" ;"Horizontal"
    resp@lbLabelStride    = 1
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@lbLabelFontHeightF   = tick_font              ; make labels smaller
    resp@lbBoxLineThicknessF  = mp_thick-2
;    resp@txFont        = font_index
;    resp@txFontHeightF = title_font 
;    resp@txFontThicknessF = title_thick
;    resp@txString      =  fig_title 
    resp@gsnPaperOrientation = "portrait"  ;"landscape"

do na = 0,narea-1,1
    res@gsnRightString  = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)" 
do nc = 0,ncase-1,1
    nn = nc + 1
    np = ncase*na + nc
    res@gsnLeftString    = pre_case(nc) + nn + ") " + case(nc)
    plot(np)  = gsn_csm_contour(wks,power(na,nc,{1:1024},:),res)
    plot2(np) = gsn_csm_contour(wks,  sig(na,nc,{1:1024},:),res_sig)
    overlay(plot(np),plot2(np))
    plot(np) = ShadeCOI(wks, plot(np), var, var&time, False)
    print(area(na))
    print(max(power(na,nc,:,:)))
    print(min(power(na,nc,:,:)))
end do
end do

gsn_panel(wks,plot,(/narea,ncase/),resp)
end

