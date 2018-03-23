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
path   = "/users/yangsong3/renql/project/TP_NUDG"
fig_name   = path + "/data/test"
filename   = path + "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.1995.daily.PS.nc" 
fig_title  = "surface pressure"
fig_out    = "pdf";"x11"
var_name   = "PS"
unit       = "Pa"

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
   mp_thick       = 5.0

;the number of the panel plot
    narea = dimsizes(area)
    plot   = new((/1/),graphic)
    line   = new((/narea/),graphic)
    ;text   = new((/narea/),graphic)

f     = addfile(filename,"r")
var1  = f->$var_name$(179,{plats:platn},{plonl:plonr}) ;get an array (nlat,nlon)
var   = var1/100
copy_VarMeta(var1,var)

;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

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
    res@mpFillOn         = False ;True
    ;res@mpLandFillColor  = "gray"
    ;res@mpFillDrawOrder  = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

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

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    ;res@cnLineThicknessF = 0
    ;res@gsnContourZeroLineThicknessF = 4.0
    ;res@cnLineColor      = "blue"
   
    res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 1          ; set min contour level
    ;res@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 65                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsEdgesOn            = True
    lnres@gsEdgeColor          = "black"
    lnres@gsEdgeDashPattern    = 0
    lnres@gsEdgeThicknessF     = 4.0
    ;fillcolor = (/"yellow","red","red","red","red","blue","blue","blue"/)
    ;opacity   = (/  1,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7 /)
    lnres@gsFillOpacityF       = 0 

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
    res@gsnRightString      = unit 
    res@gsnLeftString       = var_name ;""
    
    plot = gsn_csm_contour_map_ce(wks,var,res)
do na = 0,narea-1,1
    ;lnres@gsFillColor    = fillcolor(na)
    ;lnres@gsFillOpacityF = opacity(na)
    line(na) = gsn_add_polygon( wks, plot, lonx(:,na), latx(:,na), lnres )
    ;text(na) = gsn_add_text( wks, plot, area(na), lonm(na), latm(na),txres)
end do 

draw(plot)
frame(wks)
end

