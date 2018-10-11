;************************************************
;fourier_transform_derivative.ncl
;
;This is the fourth homework of Atmospheric Modeling
;Using discrete Fourier transform method to compute the first derivative 
;of the following functions on the interval [-1,1] with a grid spacing 2/128.
;And also output the RMS error and the plots.
;                                       by Renql
;                                      2018/10/10
;
;due to ncl do not have a complex data type, so can not calculate the fourier transform
;************************************************
begin
;************************************************
; generate variable
;************************************************
  n = 128
  f = "/home/ys17-19/renql/project/2018homework/numerical_model/derivative4.dat"
  x  = fbindirread(f,0,(/3,128/),"float")
  y  = fbindirread(f,1,(/3,128/),"float")
  dy = fbindirread(f,2,(/3,128/),"float") ;real first derivative
  ft = fbindirread(f,3,(/3,128/),"float") ;first derivative by fourier transform
  name  = (/"y=sin(pi*x)","y=exp(-5x~S~2~N~)","y=1-|x|"/)
  label = (/"y","y'","FT"/)
  
  print(x(1,10:20))
  print(y(1,10:20))
  print(dy(1,10:20))
  print(ft(1,10:20))
  printVarSummary(dy)
  if(any(ismissing(dy))) then
     print("dy contains some missing values. Beware.")
  end if
  if(any(ismissing(y))) then
     print("y contains some missing values. Beware.")
  end if
  if(any(ismissing(ft))) then
     print("ft contains some missing values. Beware.")
  end if
  plot = new((/3/),graphic)

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

wks  = gsn_open_wks("png","fourier_transform_derivative")                 ; send graphics to PNG file
  
  res               = True                        ; plot mods desired
  res@gsnFrame      = False ;True
  res@gsnDraw       = False ;True
  res@gsnMaximize   = True
  res@tiXAxisString = "x"  ; xaxis
  res@tiYAxisAngleF = 0.0 ;let the yaxis horizontal
  res@tmXBTickSpacingF   = 0.2
  res@trXMinF = -1.0
  res@trXMaxF = 1.0
  
  res@xyDashPatterns    = (/0,0,0/)
  res@xyLineThicknesses = (/4.0,4.0,4.0/)   ; make second line thicker
  res@xyLineColors      = (/"black","red","blue"/)

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
  res@lgLabelFont            = font_index
  res@pmLegendWidthF   = 0.08
  res@pmLegendHeightF  = 0.1
  res@pmLegendZone  = 0   ;set the positional origin in the center of the plot
  res@pmLegendOrthogonalPosF  = 0.4   ; move ref vector along Y-axis
  res@pmLegendParallelPosF    = -0.2    ; move ref vector along X-axis

  res@gsnStringFont         = font_index
  res@gsnStringFontHeightF = subtitle_font

;  res@tiMainString  = "Function Value"
;  res@tiYAxisString = "y"  ; xaxis
;  plot(0) = gsn_csm_xy(wks,x(0,:),y,res)   ; create plot  
  
;  res@tiMainString  = "First Derivative"
;  res@tiYAxisString = "dy"  ; xaxis
;  plot(1) = gsn_csm_xy(wks,x(0,:),dy,res)   ; create plot  

;  res@pmLegendDisplayMode = "Always"
;  res@tiMainString  = "Fourier First Derivative"
;  res@tiYAxisString = "ft"  ; xaxis
;  plot(2) = gsn_csm_xy(wks,x(0,:),ft,res)   ; create plot  
 
do ni = 0,2,1
   res@tiMainString  = name(ni) 
   plot(ni) = gsn_csm_xy(wks,x(ni,:),(/y(ni,:),dy(ni,:),ft(ni,:)/),res)   ; create plot  
end do

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"
    resp@gsnPanelRowSpec  = True

gsn_panel(wks,plot,(/1,2/),resp)
end

