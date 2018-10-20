;----------------------------------------------------------------------
;correlation_Mento_Carlo_test.ncl   2018/10/20
;----------------------------------------------------------------------
; Concepts illustrated:
;   - creat two random series ((ntest, nsample)
;   - calc the correlation between these two series 
;   - data classification to calc the probability density function
;   - draw the box plot and curve plot for cumulative probability
;----------------------------------------------------------------------
begin
;************************************************
; generate variable
;************************************************
    av = 0.0  ;representing mean of the normal distribution
    sd = 1.0  ;representing standard deviation of the normal distribution
    ntest   = 1000
    nsample = 30
    nintval = 21
    r_refer = fspan(0, 100, nintval) ;the x-axis
    frequency = new((/nintval/),float)
    cumu_freq = new((/nintval/),float)

;Generates random numbers using a normal distribution.
    x1 = random_normal(av,sd,(/ntest, nsample/))
    x2 = random_normal(av,sd,(/ntest, nsample/))
    r =  esccr(x1, x2, 0)*100 ;calc cross-corr on the rightmost dimension, get (ntest)
    r := ndtooned(r)
    printVarSummary(r)
    print(x1(1,0:10))
    print(x2(1,0:10))
    
frequency(0) = dimsizes(ind(r.eq.0))
cumu_freq(0) = frequency(0)
do ni = 1, nintval-1, 1
    frequency(ni) = dimsizes(ind((r.gt.r_refer(ni-1)).and.(r.le.r_refer(ni))))
    cumu_freq(ni) = cumu_freq(ni-1) + frequency(ni)
end do

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

wks  = gsn_open_wks("png","amplitude_error")                 ; send graphics to PNG file
  
  res               = True                        ; plot mods desired
  res@gsnFrame      = False
  res@gsnDraw       = False
  res@gsnMaximize   = True
  res@tiMainString  = "Amplitude Error"     ; title
  res@tiXAxisString = "k~F19~["  ; xaxis
  res@tiYAxisString = "A~B~k'~N~/k"; yaxis
  res@tiYAxisAngleF = 0.0 ;let the yaxis horizontal
  
  res@xyDashPatterns    = (/0,1,2/)
  res@xyLineThicknesses = (/4.0,4.0,4.0/)   ; make second line thicker
  res@xyLineColors      = (/"black","black","black"/)

  res@tmBorderThicknessF  = mp_thick
  res@tmXBLabelFont         = font_index
  res@tmYLLabelFont         = font_index
  res@tmXBLabelFontHeightF     = tick_font 
  res@tmXBLabelFontThicknessF  = tick_thick
  res@tmYLLabelFontHeightF     = tick_font
  res@tmYLLabelFontThicknessF  = tick_thick
  
  res@tmXBMajorThicknessF     = mp_thick
  res@tmYLMajorThicknessF     = mp_thick
  res@tmXBMinorThicknessF     = mp_thick
  res@tmYLMinorThicknessF     = mp_thick
  
  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 2.0           ; 2.0 is default
  res@tmXMajorGridLineColor       = "gray"
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 2.0           ; 2.0 is default
  res@tmYMajorGridLineColor       = "gray"
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  
  res@pmLegendDisplayMode = "Always"
  res@lgOrientation       = "Vertical"
  res@lgPerimOn           = False
  res@lgPerimFill         = "SolidFill"
  res@lgPerimFillColor    = "white"
  res@xyExplicitLegendLabels = label
  res@lgLabelFontColor       = "black"
  res@lgLabelFontThicknessF  = label_thick
  res@lgLabelFontHeightF     = label_font
  res@pmLegendWidthF   = 0.08
  res@pmLegendHeightF  = 0.1
  res@pmLegendZone  = 0   ;set the positional origin in the center of the plot
  res@pmLegendOrthogonalPosF  = 0.3   ; move ref vector along Y-axis
  res@pmLegendParallelPosF    = -0.35    ; move ref vector along X-axis

  res@gsnStringFont         = font_index
  res@gsnStringFontHeightF = subtitle_font

  plot = gsn_csm_xy(wks,kx,y,res)   ; create plot  

  gsres                   = True                        ; poly res
  gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
  gsres@gsFillColor       = "grey"                 ; color chosen
  xp = (/0.0,3.0,3.0,0.0,0.0/)
  yp = (/0.9,0.9,1.0,1.0,0.9/)
  dummy = gsn_add_polygon (wks,plot,xp,yp,gsres)

  draw(plot)
  frame(wks)
end


