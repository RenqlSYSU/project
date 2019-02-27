;************************************************
;shortwave_radiation.ncl
;
;This is the first homework of Tropical Dynamics
;to calc the daily solar shortwave radiation reaching the top on the atmosphere.
;The result is short_radi (lat,time), unit is J/(m2*day)
;
;                                       by Renql
;                                      2019/02/27
;************************************************

begin
;---Define parameters------------------
lat  = fspan(-90.0, 90.0, 19)
lat@long_name = "latitude"
lat@units     = "degrees_north"
print(lat(17:18))

time = ispan(1, 365, 1) ;Summer Solstice June 22 is the 173 day
nlat = dimsizes(lat)
ntime= dimsizes(time)
hour_angle = new((/nlat,ntime/),float)
short_radi = new((/nlat,ntime/),float)

pai  = 3.1415926
obliquity   = 23.45 ;unit is the degree, obliquity of the ecliptic
day_period  = 24*60*60  ;unit is s, the time of an all day
solar_const = 1365 ;W/cm2

;---calculate--------------------------------------------
day_angle   = 2*pai*(time-1)/ntime
dist_ratio  = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + \
                0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle)
sun_point   = obliquity*cos(2*pai*(173-time)/ntime) ;unit is the degree, Point of direct sunlight
print(dist_ratio(168:178))
print(dist_ratio(350:360))

do ilat = 0, nlat-1,1
    hour_angle(ilat,:) = acos(-tan(pai*sun_point/180.0)*tan(pai*lat(ilat)/180.0))
    ;short_radi(ilat,:) = (day_period/pai/10000/4.18)*solar_const*dist_ratio*(hour_angle(ilat,:)*sin(pai*lat(ilat)/180.0)*sin(pai*sun_point/180.0)+ \
    ;                     cos(pai*lat(ilat)/180.0)*cos(pai*sun_point/180.0)*sin(hour_angle(ilat,:)))
    short_radi(ilat,:) = solar_const*dist_ratio*(sin(pai*lat(ilat)/180.0)*sin(pai*sun_point/180.0)+ \
                            cos(pai*lat(ilat)/180.0)*cos(pai*sun_point/180.0))
end do
print(hour_angle(18,168:178));*180/pai)
short_radi!0    = "lat"
short_radi&lat  = lat
short_radi!1    = "time"
short_radi&time = time
unit = "cal/cm~S~2~N~/day"
;unit = "W/m~S~2~N~"
printVarSummary(hour_angle)
printVarSummary(short_radi)

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

figname = "shortwave_radiation" 
fig_title = "Incoming Solar Shortwave Radiation" 
months  = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
wks  = gsn_open_wks("png",figname)      ; send graphics to PNG file
  
    gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    res@vpHeightF        = 0.6
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiXAxisString   = "Time"
    res@tiYAxisString   = "Latitude(N)"
    
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = ispan(15,365,30)
    res@tmXBLabels      = months

    ;res@tmYLMode        = "Manual"
    ;res@tmYLTickStartF  = lats(na)
    ;res@tmYLTickEndF    = latn(na)
    ;res@tmYLTickSpacingF= lat_spacing

    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont         = font_index
    res@tmYLLabelFont         = font_index
    res@tmXBLabelFontHeightF     = tick_font 
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    res@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    res@tmXTOn               = False ; Turn off x-axis (top) major ticks
    res@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    res@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    res@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    
    res@tmXMajorGrid                   = True
    res@tmXMajorGridLineDashPattern    = 2
    res@tmXMajorGridLineColor          = "black"
    res@tmXMajorGridThicknessF         = 2.0
    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "black"
    res@tmYMajorGridThicknessF         = 2.0
    
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    ;res@cnlevelselectionmode = "automaticlevels"   ; set manual contour levels
    ;res@cnlevelselectionmode = "explicitlevels"   ; set manual contour levels
    ;res@cnlevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        =  0          ; set min contour level
    res@cnMaxLevelValF        = 750  ;32        ; set max contour level
    res@cnLevelSpacingF       =  50        ; set contour interval

    res@gsnSpreadColors      = False            ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21               ; end at color 200
    
    res@lbLabelBarOn         = True  ; the label bar
    res@lbOrientation        = "Horizontal" ;"vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 2
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-2

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;-----------set the attribute of the polyline---------------
    lres                = True
    lres@gsLineThicknessF = 7.0            ; line thickness   
    lres@gsLineColor      = "black"         ; color of lines
    
;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    
    res@tiMainString         = fig_title 
    res@tiMainFontColor      = "black"

    res@gsnLeftString   = ""
    res@gsnRightString  = unit 
    plot   = gsn_csm_lat_time(wks,short_radi,res)
    line   = gsn_add_polyline(wks,plot,time,sun_point,lres)
    print(max(short_radi))
    print(min(short_radi))

draw(plot)
frame(wks)
end

  
