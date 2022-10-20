;*********************************************************************
;filename:   draw_gtopo_map.ncl_
;input_data: gtopo30_0.9x1.25.nc 
;Draw the topography of East Asia and then label the NUDG area and
;buffer zone
;                                                     by Ql Ren
;                                                    2018/01/11
;*********************************************************************************
begin
;the path of the file
path    = "/home/ys17-19/renql/project/TP_NUDG"
filename   = path + "/analysis/mdata/gtopo30_0.9x1.25.nc"
;filename = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
;filename = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
fig_name   = path + "/analysis/fig/topography_nudg_zone_gray"
fig_title  = "CTRL";"CTRL" ; + var_name + " in " + lev + " hPa"
fig_out    = "pdf";"eps";"x11"
;fig_out    = "eps";"pdf";"x11"
var_name   = "LANDM_COSLAT"; "PHIS"
unit       = "x100 gpm"

area = (/"buffer","NUDG"/)
lats = (/20 ,25 /)
latn = (/45 ,40 /)
lonl = (/65 ,70 /)
lonr = (/110,105/)

;set the location of the area's text
latm = ( lats + latn )/2
lonm = ( lonl + lonr )/2
;lonm(5:7) = lonm(5:7) + 10 ;modify the text locate of the three areas of china
;latm(0)   = latm(0) + 5

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
   mp_thick       = 3.0
   font_index     = 22

;the number of the panel plot
    narea = dimsizes(area)
    plot   = new((/1/),graphic)
    line   = new((/narea/),graphic)
    text   = new((/2/),graphic)

f     = addfile(filename,"r")
var1  = f->$var_name$  ;({plats:platn},{plonl:plonr})$ ;get an array (nlat,nlon)
;var   = var1/9.8/100
copy_VarMeta(var1,var)

;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    gsn_define_colormap(wks,"GMT_gray") 
    ;has 11 color,0 is white, 1 is black, 2-11

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "landscape"

;-------set the attribute of the map----------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = platn
    res@mpMinLatF        = plats
    res@mpMaxLonF        = plonr
    res@mpMinLonF        = plonl
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False;True
    res@mpLandFillColor  = "gray"
    res@mpFillDrawOrder  = "PreDraw"
    ;res@mpGeophysicalLineThicknessF = mp_thick
    ;res@mpGeophysicalLineColor      = "black";"white"

;-------set the form of the coordinate--------------
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
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True;False
    res@cnLinesOn        = False;True
    res@cnLineLabelsOn   = False;True
    res@cnInfoLabelOn    = False
    res@cnLineThicknessF = 2
    ;res@gsnContourZeroLineThicknessF = 4.0
    res@cnMonoLineColor   = False
    res@cnLineColors      = (/"blue","red","green","purple"/)
    ;res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ;res@cnLevels         = (/15,25,30,40/)
    
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 10         ; set min contour level
    ;res@cnMaxLevelValF        = 50 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 5        ; set contour interval
    ;res@gsnSpreadColors      = True              ; turn of the full color map
    ;res@gsnSpreadColorStart  = 11                ; start at color 17
    ;res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsEdgesOn            = True
    lnres@gsEdgeColor          = "black"
    lnres@gsEdgeDashPattern    = 0
    lnres@gsEdgeThicknessF     = mp_thick + 1
    ;fillcolor = (/"yellow","red","red","red","red","blue","blue","blue"/)
    ;opacity   = (/  1,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7 /)
    lnres@gsFillOpacityF       = 0 

;--------set the attribute of the text------------
    txres                      = True
    txres@txFontHeightF        = subtitle_font
    txres@txFontThicknessF     = subtitle_thick
    txres@txFontColor          = "black"
    txres@txFont          = font_index 

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@tiMainString         = fig_title 
    res@tiMainFontHeightF    = title_font
    res@tiMainFontColor      = "black"
    res@gsnRightString      = unit 
    res@gsnLeftString       = var_name ;""
    
    plot = gsn_csm_contour_map_ce(wks,var,res)
    ;plot = gsn_csm_map_ce(wks,res)
do na = 0,narea-1,1
    ;lnres@gsFillColor    = fillcolor(na)
    ;lnres@gsFillOpacityF = opacity(na)
    line(na) = gsn_add_polygon( wks, plot, lonx(:,na), latx(:,na), lnres )
end do 
    ;text = gsn_add_text( wks, plot, "NUDG", lonm(0), latm(0),txres)

draw(plot)
frame(wks)
end

