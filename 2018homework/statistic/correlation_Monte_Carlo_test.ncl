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
    nsample = 10
    nintval = 21
    r_refer = fspan(0, 100, nintval) ;the x-axis
    frequency = new((/nintval/),float)
    cumu_freq = new((/nintval/),float)
    yaxis_name = (/"Cumulative Frequency","Relative Frequency"/) + "(%)"
    xaxis_name = "abs(r)*100"
    title      = "Monte_carlo_corr_test(10)"

;Generates random numbers using a normal distribution.
    x1 = random_normal(av,sd,(/ntest, nsample/))
    x2 = random_normal(av,sd,(/ntest, nsample/))
    r =  esccr(x1, x2, 0)*100 ;calc cross-corr on the rightmost dimension, get (ntest)
    r := ndtooned(abs(r))
    printVarSummary(r)
    print(x1(1,0:5))
    print(x2(1,0:5))
    
frequency(0) = dimsizes(ind(r.eq.0))
cumu_freq(0) = frequency(0)
do ni = 1, nintval-1, 1
    frequency(ni) = dimsizes(ind((r.gt.r_refer(ni-1)).and.(r.le.r_refer(ni))))
    cumu_freq(ni) = cumu_freq(ni-1) + frequency(ni)
end do
frequency = frequency/ntest*100
cumu_freq = cumu_freq/ntest*100
print(max(r))
print(min(r))
print(sum(frequency))

if(any(ismissing(r))) then
   print("r contains some missing values. Beware.")
end if

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

wks  = gsn_open_wks("png",title)                 ; send graphics to PNG file
  
  resl               = True                        ; plot mods desired
  resl@gsnFrame      = False
  resl@gsnDraw       = False
  resl@gsnMaximize   = True
  resl@vpWidthF      = 0.9
  resl@vpHeightF     = 0.3

  resl@tiMainString  = title
  resl@tiXAxisString = xaxis_name  ; xaxis
  resl@tiYAxisString = yaxis_name(0) 
  resl@tiYAxisFontColor = "black"
  
  resl@xyDashPatterns    = (/0/)
  resl@xyLineThicknesses = (/4.0/)   ; make second line thicker
  resl@xyLineColors      = (/"black"/)
  resl@trYMaxF  = 100
  resl@trYMinF  = 0
  resl@trXMaxF  = 70
  resl@trXMinF  = 0

  resl@tmBorderThicknessF  = mp_thick
  resl@tmXBLabelFont         = font_index
  resl@tmYLLabelFont         = font_index
  resl@tmXBLabelFontHeightF     = tick_font 
  resl@tmXBLabelFontThicknessF  = tick_thick
  resl@tmYLLabelFontHeightF     = tick_font
  resl@tmYLLabelFontThicknessF  = tick_thick
  
  resl@tmXBMajorThicknessF     = mp_thick
  resl@tmYLMajorThicknessF     = mp_thick
  resl@tmXBMinorThicknessF     = mp_thick
  resl@tmYLMinorThicknessF     = mp_thick
  
  resl@tmXMajorGrid                = True          ; implement x grid 
  resl@tmXMajorGridThicknessF      = 2.0           ; 2.0 is default
  resl@tmXMajorGridLineColor       = "gray"
  resl@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
  resl@tmYMajorGrid                = True          ; implement y grid
  resl@tmYMajorGridThicknessF      = 2.0           ; 2.0 is default
  resl@tmYMajorGridLineColor       = "gray"
  resl@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  
  ;resl@pmLegendDisplayMode = "Always"
  ;resl@lgOrientation       = "Vertical"
  ;resl@lgPerimOn           = False
  ;resl@lgPerimFill         = "SolidFill"
  ;resl@lgPerimFillColor    = "white"
  ;resl@xyExplicitLegendLabels = label
  ;resl@lgLabelFontColor       = "black"
  ;resl@lgLabelFontThicknessF  = label_thick
  ;resl@lgLabelFontHeightF     = label_font
  ;resl@pmLegendWidthF   = 0.08
  ;resl@pmLegendHeightF  = 0.1
  ;resl@pmLegendZone  = 0   ;set the positional origin in the center of the plot
  ;resl@pmLegendOrthogonalPosF  = 0.3   ; move ref vector along Y-axis
  ;resl@pmLegendParallelPosF    = -0.35    ; move ref vector along X-axis

  resl@gsnStringFont        = font_index
  resl@gsnStringFontHeightF = subtitle_font

  resr               = True
  resr@tiYAxisString = yaxis_name(1)
  resr@tiYAxisFontColor = "red"
  resr@xyDashPatterns    = (/0/)
  resr@xyLineThicknesses = (/4.0/)   ; make second line thicker
  resr@xyLineColors      = (/"red"/)
  plot = gsn_csm_xy2(wks,r_refer,cumu_freq,frequency,resl,resr)   ; create plot  

  gsres                   = True                        ; poly res
  gsres@tfPolyDrawOrder  = "Predraw"                   ; draw this first
  gsres@gsFillColor      = "grey"                 ; color chosen
  xp = (/0.0,100,100,0.0,0.0/)
  yp = (/90 ,90 ,100,100,90 /)
  dummy = gsn_add_polygon (wks,plot,xp,yp,gsres)

  draw(plot)
  frame(wks)
end


