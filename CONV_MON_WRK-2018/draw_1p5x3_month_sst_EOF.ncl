;******************************************************************************
;filename:      draw_1p5x3_month_sst_EOF.ncl
;input_file:    var_name(nv) + area + _EOFn3.nc
;----------------------------------------------------------------------
; Concepts illustrated:
;   - read EOF pattern and time series data
;   - Draw the first three EOF spatial patterns and principal components
;   - each panel have four map figure and one time series figure(have four curves)
;                                            by Ql_Ren
;                                           2018/09/17
;*********************************************************************************
begin
unit     = (/"~S~o~N~C","hPa","W/m2","W/m2"/)    ;this can be used as right string 
var_name = (/"TS","PSL","FLNT","FSNT"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)
nv = 0

;-----EOF area-------------- 
lats = -30
latn = 30
lonl = 90
lonr = 300
;lev  = 200
neval = 4
lat_spacing = 15
lon_spacing = 20

opt = True
opt@jopt = 0 ;0 mean use the covariance matrix to compute EOFs (default), 1 mean use the correlation matrix
eof_kind = (/"cov","corr"/)
tname    = (/"_","_DJF_"/)

;----path of the file and output file-----------------------------
path      = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename  = path + "/mdata/" + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E" + tname(0) + eof_kind(opt@jopt) + "EOF" + neval + ".nc"
fig_name  = path + "/fig/"   + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E" + tname(0) + eof_kind(opt@jopt) + "EOF" + neval
fig_out   = "pdf";"x11"

;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
map_var = f->eof       ;(ncase,neval,nlat,nlon)
ts_var  = dim_standardize_n_Wrap(f->eof_ts,1,2)    ;(ncase,neval,nyear) 
right_title   = f->pcvar
center_title  = f->eval
printVarSummary(map_var)
printVarSummary(ts_var)
printVarSummary(right_title)

ncase = dimsizes(ts_var&case)
ntime = dimsizes(ts_var&time)
sig   = new((/ncase,neval/),logical) 
plot  = new((/(ncase+1)*(neval-1)/),graphic)

do nc = 0,ncase-1,1
sig(nc,:) = eofunc_north(right_title(nc,:),ntime,True)
end do
subtitle_color = where(sig,"red","black")

;==============================================================
;draw figure 
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.02
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"

;------------set the attribute of the map----------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = lat_spacing 
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;------------set the form of the coordinate-----------------------
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

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -0.063       ; set min contour level
    res@cnMaxLevelValF        = 0.063       ; set max contour level
    res@cnLevelSpacingF       = 0.007        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    
;------------set the resxy for time series figure-----------
    resxy     = True
    ;resxy@vpXF      = 0.1
    ;resxy@vpYF      = 0.4
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    
    resxy@xyDashPatterns     = (/0,0,0,0/)
    resxy@xyLineThicknesses = (/2.0, 2.0, 2.0, 2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black", "red", "green", "blue"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    resxy@trYMinF = -3.0
    resxy@trYMaxF = 3.0
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = label_font
    resxy@tmXBLabelFontThicknessF  = label_thick 
    resxy@tmYRLabelFontHeightF     = label_font
    resxy@tmYRLabelFontThicknessF  = label_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    resxy@tmXBMode           = "Manual"
    resxy@tmYLMode           = "Manual"
    resxy@tmXBTickSpacingF   = ntime/10
    ;resxy@tmYLTickSpacingF   = 1.0

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    resxy@gsnYRefLine             = 0
    resxy@gsnYRefLineDashPatterns = 0 ; (/12,12/)   
    resxy@gsnYRefLineThicknesses  = 3.0 ;(/3.0,3.0/)
    resxy@gsnYRefLineColors       = "black"  ;(/"blue","red"/)
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Horizontal"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case 
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFont            = font_index
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.8
    resxy@pmLegendHeightF  = 0.03
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = -0.6   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0    ; move ref vector along X-axis

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txFont = font_index
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    ;resp@gsnPanelRowSpec  = True
    ;resp@gsnPanelBottom   = 0.2
    resp@gsnPaperOrientation = "landscape" ;"portrait"

do nc = 0,ncase-1,1
do ne = 0,neval-2,1
    np = nc*(neval-1) + ne
    res@gsnLeftString    = case(nc)
    res@gsnCenterString  = "EOF " + (ne+1)
    res@gsnRightStringFontColor    = subtitle_color(nc,ne)
    res@gsnRightString   = center_title(nc,ne)+"  "+right_title(nc,ne)+"%"
    if(nc.le.1) then
        res@mpGeophysicalLineThicknessF = 0
    else
        res@mpGeophysicalLineThicknessF = mp_thick
    end if
    plot(np) = gsn_csm_contour_map_ce(wks,map_var(nc,ne,:,:),res)
    print(max(map_var(nc,ne,{lats:latn},{lonl:lonr})))
    print(min(map_var(nc,ne,{lats:latn},{lonl:lonr})))
end do
end do 
    
do ne = 0,neval-2,1
    np = 12 + ne
    plot(np) = gsn_csm_xy(wks,ts_var&time,(/ts_var(:,ne,:)/),resxy)
end do 

    resp@txString      = "EOF of " + var_name(nv)
    gsn_panel(wks,plot,(/5,3/),resp)
end

