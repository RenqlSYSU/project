;*********************************************************************
;filename:      draw_1p2x1_preci_wind_obs_ctrl.ncl
;use cpc global precip daily 1x1 data and NCEP2 analysis wind data 
;to state the CTRL's results are right
;                                                     by Ql Ren
;                                                    2018/0418
;*********************************************************************************
begin
draw_var = (/"500hPa Wind & Preci(mm/day)"/)  ;this is the variable used to draw, can used as center string
;draw_varname = (/"preci","uwnd","vwnd","preci","U","V"/)
draw_varname = (/"preci","uwnd","vwnd","preci","U","V"/)
var_name = (/"season_ave"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
lev      = (/850,500,200/)
nm = 2
nl = 1;0

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/analysis/mdata/" + (/"Observ","CTRL"/) + "-Clim_month_ave_" ;+ draw_varname(nv) + ".nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/analysis/fig/1p2x1_preci_wind_obs_ctrl_" + months(nm) + lev(nl)
    case   = (/" Obs"," CTRL"/)
pre_case   = (/"(a)","(b)"/)
fig_out    = "eps";"pdf";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8
land   = f->LANDM_COSLAT({lats:latn},{lonl:lonr})
ncase  = dimsizes(case)

;the number of the panel plot
    plot     = new((/ncase/),graphic)
    line     = new((/ncase/),graphic)
    gtop     = new((/ncase/),graphic)
    plot_uv   = new((/ncase/),graphic)

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
   vcref_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 3.0
   mp_thick       = 5.0
   font_index     = 22

;================================================================
;read data
;====================================================================
nf = 0
f = addfile(filename(nf)+draw_varname(0)+".nc","r")
 var0  = f->$draw_varname(0)$({lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
 ;var0  = f->$draw_varname(0)$(nl,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
f = addfile(filename(nf)+draw_varname(1)+".nc","r")
 var0_u  = f->$draw_varname(1)$(nl,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
f = addfile(filename(nf)+draw_varname(2)+".nc","r")
 var0_v  = f->$draw_varname(2)$(nl,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)

nf = 1
f = addfile(filename(nf)+draw_varname(3)+".nc","r")
 var1  = f->$var_name$(nm,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
 ;var1  = f->$var_name$(nm,nl,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
f = addfile(filename(nf)+draw_varname(4)+".nc","r")
 var1_u  = f->$var_name$(nm,nl,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
f = addfile(filename(nf)+draw_varname(5)+".nc","r")
 var1_v  = f->$var_name$(nm,nl,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)

;var1 = mask(var1,land.lt.1,False)
printVarSummary(var0)
printVarSummary(var0_u)
printVarSummary(var0_v)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 10
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
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    
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
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 0          ; set min contour level
    res@cnMaxLevelValF        = 15         ; set max contour level
    res@cnLevelSpacingF       = 1          ; set contour interval

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "black"
    resvc@vcMinDistanceF             = 0.03
    resvc@vcRefMagnitudeF            = 8.0             ; define vector ref mag
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = -0.67           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = -0.1            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    resvc@vcRefAnnoString2           = "8m/s"
    
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "
    resvc@vcLineArrowHeadMinSizeF   = 0.01
    resvc@vcLineArrowHeadMaxSizeF   = 0.01

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 4.0

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
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString  = draw_var   

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"
    resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    resp@lbLabelFont           = font_index              ; make labels smaller
    resp@lbLabelFontHeightF    = tick_font-0.005
    resp@lbLabelFontThicknessF = tick_thick
    resp@lbLabelStride         = 1
    resp@lbBoxLineThicknessF   = mp_thick-2

nc = 0
res@gsnLeftString   = pre_case(nc) + case(nc)
plot(nc)     = gsn_csm_contour_map_ce(wks,var0,res)
plot_uv(nc)  = gsn_csm_vector(wks,var0_u,var0_v,resvc)
overlay(plot(nc),plot_uv(nc))
gtop(nc) = gsn_csm_contour(wks,phis,res_gtop)
overlay(plot(nc),gtop(nc))
line(nc) = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
print(max(var0))
print(min(var0))

nc = 1
res@gsnLeftString   = pre_case(nc) + case(nc)
plot(nc)     = gsn_csm_contour_map_ce(wks,var1,res)
plot_uv(nc)  = gsn_csm_vector(wks,var1_u,var1_v,resvc)
overlay(plot(nc),plot_uv(nc))
gtop(nc) = gsn_csm_contour(wks,phis,res_gtop)
overlay(plot(nc),gtop(nc))
line(nc) = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
print(max(var1))
print(min(var1))

gsn_panel(wks,plot,(/1,ncase/),resp)
end

