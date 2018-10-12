;*****************************************************************
;phase_group_velocity_error.ncl
;
;This is the fourth homework of Atmospheric Modeling, which asks us to 
;visualize the phase velocity error and group velocity error
;of the 2-order CTD,4-order CTD and 4-order CTFD
;                                           by Renql
;                                         2018/10/12
;****************************************************************
begin
;************************************************
; generate variable
;************************************************
  kx  = fspan(0.001,6,60)
  nkx = dimsizes(kx)
  c   = new((/3,nkx/),float)
  cg  = new((/3,nkx/),float)
  c(0,:)  = sin(kx)/kx                ;c for 2-order CTD
  cg(0,:) = cos(kx)                   ;cg for 2-order CTD
  c(1,:)  = sin(kx)*(4.0-cos(kx))/3.0/kx                 ;c for 4-order CTD
  cg(1,:) = 4.0/3*cos(kx)-1.0/3*cos(2*kx)                ;cg for 4-order CTD
  c(2,:)  = 3*sin(kx)/(2+cos(kx))/kx                   ;c for 4-order CTFD
  cg(2,:) = (3+6*cos(kx))/(2+cos(kx))^2                ;cg for 4-order CTFD

  label  = (/"c' 2-order CTD","c' 4-order CTD","c' 4-order CTFD",\
             "c~B~g~N~' 2-order CTD","c~B~g~N~' 4-order CTD","c~B~g~N~' 4-order CTFD"/)

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

wks  = gsn_open_wks("png","phase_group_velocity_error")      ; send graphics to PNG file
  
  res               = True                        ; plot mods desired
  res@gsnFrame      = False
  res@gsnDraw       = False
  res@gsnMaximize   = True
  res@tiXAxisString = "k~F19~["  ; xaxis
  res@tiYAxisAngleF = 0.0 ;let the yaxis horizontal
  
  res@xyDashPatterns    = (/0,1,2,0,1,2/)
  res@xyLineThicknesses = (/4.0,4.0,4.0,4.0,4.0,4.0/)   ; make second line thicker
  res@xyLineColors      = (/"red","red","red","blue","blue","blue"/)

  res@gsnYRefLine             = 0.0
  res@gsnYRefLineDashPatterns = 11   
  res@gsnYRefLineThicknesses  = 4.0
  res@gsnYRefLineColors       = "black"

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
  res@pmLegendHeightF  = 0.2
  res@pmLegendZone  = 0   ;set the positional origin in the center of the plot
  res@pmLegendOrthogonalPosF  = 0.3   ; move ref vector along Y-axis
  res@pmLegendParallelPosF    = -0.3    ; move ref vector along X-axis

  res@gsnStringFont         = font_index
  res@gsnStringFontHeightF = subtitle_font

  res@tiMainString  = "Velocity Error"     ; title
  res@tiYAxisString = "1/c"; yaxis
  plot = gsn_csm_xy(wks,kx,(/c(0,:),c(1,:),c(2,:),cg(0,:),cg(1,:),cg(2,:)/),res)   ; create plot  

  draw(plot)
  frame(wks)
end


