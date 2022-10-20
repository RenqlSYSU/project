;*********************************************************************
;filename:   draw_area_descrip.ncl
;input_data:  
;Draw the global map and then label the areas ("NUDG","India","Tibet",
;"Philippines","Maritime_continent","North China","Central East China","South China")
;                                                     by Ql Ren
;                                                    2018/01/11
;*********************************************************************************
begin
;the path of the file
fig_name   = "/users/yangsong3/renql/project/TP_NUDG/nudgtest/fig/area_descrip"
fig_title  = "Area Description" ; + var_name + " in " + lev + " hPa"
fig_out    ="pdf";"x11"

area  = (/"NUDG","India","Tibet","Philippines","Maritime_continent",\
          "North China","Central East China","South China"/)
lats = (/25 ,20 ,25  ,5   ,-10 ,34  ,28  ,20 /)
latn = (/40 ,35 ,35  ,20  ,7.5 ,43  ,34  ,28 /)
lonl = (/70 ,60 ,90  ,120 ,85  ,107 ,107 ,107/)
lonr = (/105,80 ,103 ,140 ,120 ,122 ,122 ,122/)

;set the location of the area's text
latm = ( lats + latn )/2
lonm = ( lonl + lonr )/2
lonm(5:7) = lonm(5:7) + 10 ;modify the text locate of the three areas of china
latm(0)   = latm(0) + 5

;set the location of the box
lonx = (/ lonl, lonr, lonr, lonl, lonl /) ;get (5,narea)
latx = (/ lats, lats, latn, latn, lats /)

;plot area of Asia
plats = -15
platn = 55
plonl = 30
plonr = 150

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 6.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0

;the number of the panel plot
    narea = dimsizes(area)
    plot   = new((/1/),graphic)
    line   = new((/narea/),graphic)
    text   = new((/narea/),graphic)

;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    ;res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "landscape"

;-------set the attribute of the map----------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = platn
    res@mpMinLatF        = plats
    res@mpMaxLonF        = plonr
    res@mpMinLonF        = plonl
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = True
    res@mpLandFillColor = "gray"
    res@mpFillDrawOrder = "PreDraw"
    ;res@mpGeophysicalLineThicknessF = mp_thick
    ;res@mpGeophysicalLineColor      = "black";"white"

;-------set the form of the coordinate--------------
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsEdgesOn            = True
    lnres@gsEdgeColor          = "black"
    lnres@gsEdgeDashPattern    = 0
    lnres@gsEdgeThicknessF     = 4.0
    fillcolor = (/"yellow","red","red","red","red","blue","blue","blue"/)
    opacity   = (/  1,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7 /)

;--------set the attribute of the text------------
    txres                      = True
    txres@txFontHeightF        = subtitle_font
    txres@txFontThicknessF     = subtitle_thick
    txres@txFontColor          = "black"

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@tiMainString         = fig_title 
    res@tiMainFontHeightF    = title_font
    res@tiMainFontColor      = "black"
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    res@gsnRightString      = "" 
    res@gsnLeftString       = ""
    
    plot = gsn_csm_map_ce(wks,res)
do na = 0,narea-1,1
    lnres@gsFillColor    = fillcolor(na)
    lnres@gsFillOpacityF = opacity(na)
    line(na) = gsn_add_polygon( wks, plot, lonx(:,na), latx(:,na), lnres )
    text(na) = gsn_add_text( wks, plot, area(na), lonm(na), latm(na),txres)
end do 

draw(plot)
frame(wks)
end

