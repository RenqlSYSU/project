;************************************************
;amplitude_error_visualize.ncl
;
;This is the third homework of Atmospheric Modeling, which asks us to 
;visualize the amplitude error and phase shift of the CTFD, FWD, CTD
;************************************************
begin
;************************************************
; generate variable
;************************************************
  kx  = fspan(0.001,3,30)
  nkx = dimsizes(kx)
  y   = new((/3,nkx/),float)
  y(0,:) = 3*sin(kx)/(2+cos(kx))/kx  ;CTFD
  y(1,:) = 2*sin(kx/2)/kx            ;FWD
  y(2,:) = sin(kx)/kx                ;CTD
  label  = (/"CTFD","FWD","CTD"/)

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


