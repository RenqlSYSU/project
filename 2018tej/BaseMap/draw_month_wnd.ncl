;******************************************************************************
;filename:    draw_month_wnd.ncl
;input_data:  clim.mon.uvwnd.nc 
;draw the climate monthly wind of U and V 
;and plot the contour of U wind to indicate tej
;                                            by Ql_Ren
;                                           2016/01/02
;**********************************************************************
begin
;the filename of the read-in-data and figure, as well as the title of the figure
filename   = "/home/ys17-19/renql/project/2018tej/BaseMap/mdata/clim.mon.uvwnd.nc"
fig_name   = "/home/ys17-19/renql/project/2018tej/BaseMap/fig/12panel_3X1_month_wnd" 
fig_title  = "Climate Monthly Wind"
left_title = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;right_title= lev+"hPa"
fig_out    = "pdf";"x11"
var_name   = (/"uwnd","vwnd"/)

;the number of the panel plot
plot   = new((/3/),graphic)

;area: east Asia
lev   = (/200,150,100/)  ;the hpa of the level
;lats = -20
;latn = 60
;lonl = 20
;lonr = 220

;define the thickness and height of the title,tick and label
tick_font      = 0.02
subtitle_font  = 0.02
title_font     = 0.025
vcref_font     = 0.01
tick_thick     = 8.0
subtitle_thick = 8.0  ;not have this attribute
title_thick    = 8.0
mp_thick       = 5.0
vc_thick       = 1.5

;================================================================
;read data
;====================================================================
f     = addfile(filename,"r")
var_u = f->$var_name(0)$(:,{lev},:,:)
var_v = f->$var_name(1)$(:,{lev},:,:)
printVarSummary(var_u)
printVarSummary(var_v)

;==============================================================
;draw figure in twelve months (one month is a panel)
;==================================================================
wks = gsn_open_wks(fig_out,fig_name)
do m = 0,11,1
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True
    res@gsnPaperOrientation = "portrait"

;--------set the attribute of the map---------------
    res@mpCenterLonF     = 180.
    ;res@mpMaxLatF        = latn
    ;res@mpMinLatF        = lats
    ;res@mpMaxLonF        = lonr
    ;res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 20
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate----------------
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

;--------set the attribute of the contour--------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF       = 0          ; set min contour level
    res@cnMaxLevelValF       = 32          ; set max contour level
    res@cnLevelSpacingF      = 4         ; set contour interval

    res@gsnSpreadColors      = True               ; turn of the full color map
    res@gsnSpreadColorStart = 2                ; start at color 17
    res@gsnSpreadColorEnd   = 11                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1

;----------set the attribute of the vector-------------------
     resv@vcGlyphStyle               = "CurlyVector"
     resv@vcLineArrowThicknessF      = vc_thick
     resv@vcMinDistanceF             = 0.015
     ;resv@vcRefMagnitudeF            = 30.0             ; define vector ref mag
     resv@vcRefLengthF               = 0.05           ; define length of vec ref
     resv@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
     ;resv@vcRefAnnoSide              = "Bottom"          
     ;resv@vcRefAnnoJust              = "TopRight"
     resv@vcRefAnnoOrthogonalPosF    = 0.5           ; move ref vector along Y-axis
     resv@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
     resv@vcRefAnnoPerimOn           = False         ;close the border
     resv@vcRefAnnoFontHeightF       = vcref_font
     resv@vcRefAnnoString1On         = False
     resv@vcRefAnnoString2On         = True
     ;resv@vcRefAnnoString2           = "30m/s"

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; resv@tiMainString         = fig_title + "in" + months(m)
    ;resv@tiMainFontHeightF    = 0.04
   ; resv@tiMainFontColor      = "black"
    resv@gsnCenterString      = ""
    resv@gsnRightString       = right_title
    resv@gsnRightStringFontHeightF = subtitle_font
    resv@gsnLeftStringFontHeightF  = subtitle_font
    
    do i = 0, nvar-1,1
        res@gsnLeftString    = left_title(i)
        ;if(i.eq.2) then
          res@gsnSpreadColors = False
          res@cnMinLevelValF        = minlevel(h,i)         ; set min contour level
          res@cnMaxLevelValF        = maxlevel(h,i)          ; set max contour level
          res@cnLevelSpacingF       = levelspacing(h,i)          ; set contour interval
          res@vcRefMagnitudeF       = vcref(h,i)             ; define vector ref mag
          res@vcRefAnnoString2      = vcref(h,i)+"m/s"
        ;end if
        ;plot(i) = gsn_csm_contour_map_ce(wks,var(i,m,:,:),res)
        plot(i) = gsn_csm_vector_scalar_map_ce(wks,var(i,0,m,:,:),\
                    var(i,1,m,:,:),var(i,2,m,:,:),res)
    print(max(var(i,2,m,{lats:latn},{lonl:lonr})))
    print(min(var(i,2,m,{lats:latn},{lonl:lonr})))
    end do
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title +" " + lev(h) + " in " + months(m)
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    gsn_panel(wks,plot,(/3,1/),resp)
;draw(plot)
;frame(wks)
end do ;the loop of the month
end do ;the loop of the levels
end

