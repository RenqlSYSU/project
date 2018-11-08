;*********************************************************************
;filename:      draw_12p2x1_TN01_wave_flux_dclim.ncl
;input_data:    CTRL_bNUDG6h_TN01_wave_flux_dclim_seasonly.nc
;***********************************************************************
;read Fx,Fy,Fz to draw vector and contour plot
;                                                     by Ql Ren
;                                                    2018/11/08
;*********************************************************************************
begin
draw_var = (/"Fx","Fy","Fz"/)  ;this is the variable used to draw, can used as center string
unit     = "m~S~2~N~/s~S~2~N~" ;this can be used as right string 
months   = (/"DJF","MAM","JJA","SON"/)
pre_case = (/"(a)","(b)","(c)","(d)"/)
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
lev  = (/850,500,200/)  ;can used as right string
nl = 1;1

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
case = (/"CTRL","NUDG24h","NUDG6h"/)
basic = 2
filename = path + "/energy-20180417/mdata/CTRL_b" + case(basic) + "_TN01_wave_flux_dclim_seasonly.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/analysis/fig/4p3X1_CTRL_b" + case(basic) + "TN01_wave_flux_dclim_seasonly_"+ lev(nl)
fig_title  = "TN01 Wave Flux"
fig_out    ="pdf";"x11"

;area: east Asia
lats = -90 ;-20
latn = 90  ;60
lonl = 0  ;20
lonr = 360 ;220
siglvl = 0.05

nmonth = dimsizes(months)
ftop = addfile(gtopfile,"r")
phis = ftop->PHIS({lats:latn},{lonl:lonr})
phis = phis/9.8   ;convert unit from m2/s2 to gpm
f = addfile(filename(0),"r")
var_u  = f->$draw_var(0)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
var_v  = f->$draw_var(1)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
var_w  = f->$draw_var(2)$(:,{lev(nl)},{lats:latn},{lonl:lonr})

;the number of the panel plot
    plot     = new((/nmonth/),graphic)
    line     = new((/nmonth/),graphic)
    gtop     = new((/nmonth/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   vcref_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   vc_thick       = 1.5
   font_index     = 22

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 20
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
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
    res@gsnScalarContour = True
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -8          ; set min contour level
    res@cnMaxLevelValF        = 8 ;32        ; set max contour level
    res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

;----------set the attribute of the vector-------------------
    res@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    res@vcVectorDrawOrder          = "PostDraw"
    res@vcLineArrowThicknessF      = vc_thick
    res@vcMonoLineArrowColor       = True
    res@vcLineArrowColor           = "black"
    res@vcMinDistanceF             = 0.02
    ;res@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    res@vcRefLengthF               = 0.05           ; define length of vec ref
    res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;res@vcRefAnnoSide              = "Bottom"          
    ;res@vcRefAnnoJust              = "TopRight"
    res@vcRefAnnoOrthogonalPosF    = 0.5           ; move ref vector along Y-axis
    res@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    res@vcRefAnnoPerimOn           = False         ;close the border
    res@vcRefAnnoFont              = font_index
    res@vcRefAnnoFontHeightF       = vcref_font
    res@vcRefAnnoString1On         = False
    res@vcRefAnnoString2On         = True
    ;res@vcRefAnnoString2           = "30m/s"

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 4.0
    ;lnres@gsLineLabelString    = "TP"

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = 4.0
    res_gtop@cnLineColor      = "black"
    res_gtop@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_gtop@cnLevels         = (/1500/)
    res_gtop@gsnCenterString      = ""
    res_gtop@gsnLeftString    = " "
    res_gtop@gsnRightString   = " "

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    vcref         = (/(/100,100,100,100/), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/100,100,100,100/), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/15,15,15,15/)/)   ;the level of 200hpa    ctrl,nudg,diff
do nm = 0 ,nmonth-1,1
    res@gsnLeftString   = pre_case(nm) + months(nm)
    res@vcRefMagnitudeF       = vcref(nl,nm)                ; define vector ref mag
    res@vcRefAnnoString2      = vcref(nl,nm)+"m/s"
    plot(nm) = gsn_csm_vector_scalar_map_ce(wks,var_u(nm,:,:),var_v(nm,:,:),var_w(nm,:,:),res)
    gtop(nm) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(nm),gtop(nm))
    line(nm) = gsn_add_polyline( wks, plot(nm), lonx, latx, lnres )
    print(max(var_u(nm,{lats:latn},{lonl:lonr})))
    print(min(var_u(nm,{lats:latn},{lonl:lonr})))
    print(max(var_w(nm,{lats:latn},{lonl:lonr})))
    print(min(var_w(nm,{lats:latn},{lonl:lonr})))
end do 

    resp@txString      = fig_title
    gsn_panel(wks,plot,(/2,2/),resp)
end

