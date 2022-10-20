;*********************************************************************************
;filename:      draw_3panel_4X1_JJAu_EOF.ncl
;input_data:    Asia_JJAu_covEOFn3.nc
;Draw the first three EOF spatial patterns and principal components of Asia JJA preci data
;every panel have three figures of spatial patterns (CTRL,NUDG,DIFF) 
;and one figure of time series which has two curves (CTRL,NUDG).
;                                            by Ql_Ren
;                                           2017/12/19
;*********************************************************************************
begin
;the path of the file
filename   = "/users/yangsong3/Model_Team/F/clm-month/mdata/Asia_JJAv_covEOFn3.nc"
fig_name   = "/users/yangsong3/Model_Team/F/clm-month/fig/3panel_4X1_JJAv_EOF"
fig_title  = "of Asia JJA V(m/s) at 200 hPa"  ;"EOF" + neval + fig_title
left_title = (/"AMIP_CTRL","TP_NUDG","AMIP_CTRL - TP_NUDG"/)
fig_out    = "pdf";"x11"

;plot area 
lats = 30
latn = 50
lonl = 30
lonr = 130

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.02
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 5.0

;the number of the panel plot
    plot   = new((/4/),graphic)
    sig    = new((/2,4/),logical) 
;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
map_var = f->eof       ;(ncase,neval,nlat,nlon)
ts_var  = f->eof_ts    ;(ncase,neval,nyear) 
right_title   = f->pcvar
center_title  = f->eval
neval   = 3
printVarSummary(map_var)
printVarSummary(ts_var)
printVarSummary(right_title)
diff   = map_var(0,:,:,:) - map_var(1,:,:,:)
copy_VarMeta(map_var(0,:,:,:),diff)

do nc = 0,1
sig(nc,:) = eofunc_north(right_title(nc,:),27,True)
end do
subtitle_color = where(sig,"red","black")

;==============================================================
;draw figure 
;==================================================================
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
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 10
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;------------set the form of the coordinate-----------------------
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

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -0.09       ; set min contour level
    res@cnMaxLevelValF        = 0.09       ; set max contour level
    res@cnLevelSpacingF       = 0.01        ; set contour interval

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
    ;resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns     = (/0,0/)
    resxy@xyLineThicknesses = (/2.0, 2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"blue","red"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    resxy@trYMinF = -3.0
    resxy@trYMaxF = 3.0
    
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
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    resxy@tmXBMode           = "Manual"
    resxy@tmYLMode           = "Manual"
    resxy@tmXBTickSpacingF   = 3
    ;resxy@tmYLTickSpacingF   = 1.0

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    resxy@gsnYRefLine             = 0
    resxy@gsnYRefLineDashPatterns = 0 ; (/12,12/)   
    resxy@gsnYRefLineThicknesses  = 3.0 ;(/3.0,3.0/)
    resxy@gsnYRefLineColors       = "black"  ;(/"blue","red"/)
                             
    resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
    resxy@gsnRightStringFontHeightF  = subtitle_font
    resxy@gsnLeftStringFontHeightF   = subtitle_font
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = left_title(0:1)
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.1
    resxy@pmLegendHeightF  = 0.05
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.38   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.4    ; move ref vector along X-axis

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    ;res@gsnCenterStringFontHeightF= subtitle_font

do ne = 0,2,1
    do nc = 0,2,1
        res@gsnLeftString    = left_title(nc)
        if(nc.eq.2) then
        res@gsnRightString   = "" 
        plot(nc) = gsn_csm_contour_map_ce(wks,diff(ne,:,:),res)
        print(max(diff(ne,{lats:latn},{lonl:lonr})))
        print(min(diff(ne,{lats:latn},{lonl:lonr})))
        else
        res@gsnRightStringFontColor    = subtitle_color(nc,ne)
        res@gsnRightString   = center_title(nc,ne)+"  "+right_title(nc,ne)+"%" 
        plot(nc) = gsn_csm_contour_map_ce(wks,map_var(nc,ne,:,:),res)
        print(max(map_var(nc,ne,{lats:latn},{lonl:lonr})))
        print(min(map_var(nc,ne,{lats:latn},{lonl:lonr})))
        end if
    end do
    
    resxy@gsnCenterString  = ""
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = "" 
    ;resxy@gsnYRefLine      = (/ave(nv,0,ns,nl),ave(nv,1,ns,nl)/)
    plot(3) = gsn_csm_xy(wks,ts_var&year,(/ts_var(:,ne,:)/),resxy)

    resp  = True    
    resp@gsnMaximize   = True
    ineval=ne+1
    resp@txString      = "EOF" + ineval + " " + fig_title 
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    ;resp@gsnPanelRowSpec  = True
    ;resp@gsnPanelBottom   = 0.2
    resp@gsnPaperOrientation = "portrait"
    gsn_panel(wks,plot,(/4,1/),resp)
end do 
end

