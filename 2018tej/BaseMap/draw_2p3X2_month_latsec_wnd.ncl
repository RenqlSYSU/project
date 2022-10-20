;******************************************************************************
;filename:    draw_2p3X2_month_latsec_wnd.ncl
;input_data:  clim.mon.uvwnd.nc, clim.mon.omega.nc 
;draw the climate monthly wind of vwnd and omega
;and plot the contour of uwnd to indicate tej
;                                            by Ql_Ren
;                                           2016/01/02
;**********************************************************************
begin
;the filename of the read-in-data and figure, as well as the title of the figure
filename   = "/home/ys17-19/renql/project/2018tej/BaseMap/mdata/clim.mon.uvwnd.nc"
filename2  = "/home/ys17-19/renql/project/2018tej/BaseMap/mdata/clim.mon.omega.nc"
fig_name   = "/home/ys17-19/renql/project/2018tej/BaseMap/fig/2p3X2_month_latsec_wnd" 
fig_title  = "Climate Monthly Wind along 0-160E"
left_title = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;right_title= lev+"hPa"
fig_out    = "pdf";"x11"
var_name   = (/"uwnd","vwnd","omega"/)

;the number of the panel plot
plot   = new((/12/),graphic)
line   = new((/12/),graphic)

;area: east Asia
;lev   = (/150/)  ;the hpa of the level
;lats = 0
;latn = 160
lonl = 0
lonr = 160

;define the thickness and height of the title,tick and label
tick_font      = 0.03
subtitle_font  = 0.015
title_font     = 0.025
vcref_font     = 0.02
tick_thick     = 8.0
subtitle_thick = 8.0  ;not have this attribute
title_thick    = 8.0
mp_thick       = 4.0
vc_thick       = 2.5

;================================================================
;read data
;====================================================================
f     = addfile(filename,"r")
var_u = f->$var_name(0)$(:,{1000:100},:,{lonl:lonr})
var_v = f->$var_name(1)$(:,{1000:100},:,{lonl:lonr})
f     = addfile(filename2,"r")
var_w = f->$var_name(2)$(:,{1000:100},:,{lonl:lonr})

var_u := dim_avg_n_Wrap(var_u,3)
var_v := dim_avg_n_Wrap(var_v,3)
var_w := dim_avg_n_Wrap(var_w,3)
var_w  = -var_w*100
printVarSummary(var_u)
printVarSummary(var_v)
printVarSummary(var_w)
levx = var_u&level
latx = new(dimsizes(levx),float)
latx = 0.0

;==============================================================
;draw figure in twelve months (one month is a panel)
;==================================================================
wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnPaperOrientation = "portrait"
    res@gsnLeftString       = ""
    res@gsnRightString      = ""
    res@gsnCenterString      = ""
    
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
    res@tmYRMode                 = "Automatic"
    res@tiYAxisString            = ""

;--------set the attribute of the contour--------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
    res@cnInfoLabelOn    = False
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF       = -24          ; set min contour level
    res@cnMaxLevelValF       = 24         ; set max contour level
    res@cnLevelSpacingF      = 3         ; set contour interval

    res@gsnSpreadColors     = False ;True               ; turn of the full color map
    res@gsnSpreadColorStart = 2                ; start at color 17
    res@gsnSpreadColorEnd   = 11                ; end at color 200

    res@lbLabelBarOn         = False ;True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1

;----------set the attribute of the vector-------------------
     res@vcGlyphStyle               = "CurlyVector"
     res@vcLineArrowThicknessF      = vc_thick
     res@vcLineArrowHeadMaxSizeF    = 0.01
     res@vcLineArrowHeadMinSizeF    = 0.01
     res@vcMinDistanceF             = 0.015
     res@vcRefMagnitudeF            = 3.0             ; define vector ref mag
     res@vcRefLengthF               = 0.05           ; define length of vec ref
     res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
     ;res@vcRefAnnoSide              = "Bottom"          
     ;res@vcRefAnnoJust              = "TopRight"
     res@vcRefAnnoOrthogonalPosF    = 0.4           ; move ref vector along Y-axis
     res@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
     res@vcRefAnnoPerimOn           = False         ;close the border
     res@vcRefAnnoFontHeightF       = vcref_font
     res@vcRefAnnoString1On         = False
     res@vcRefAnnoString2On         = True
     res@vcRefAnnoString2           = "3m/s"

;---------set the attribute of the line----------------
    plres                   = True
    plres@gsLineColor       = "black"
    plres@gsLineThicknessF  = mp_thick

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; resv@tiMainString         = fig_title + "in" + months(m)
    ;resv@tiMainFontHeightF    = 0.04
   ; resv@tiMainFontColor      = "black"
    ;res@gsnRightString       = right_title
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    
do np = 0,1,1
    nm1 = 6*np
    nm2 = 6*np+5
    do nm = nm1, nm2, 1
        plot(nm) = gsn_csm_pres_hgt_vector(wks,var_u(nm,{1000:100},{-90:90}),var_v(nm,{1000:100},{-90:90}),var_w(nm,{1000:100},{-90:90}),res)
        line(nm) = gsn_add_polyline(wks,plot(nm),latx,levx,plres)
        print(min(var_u(nm,:,:)))
    end do
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    resp@gsnPanelFigureStrings = left_title(nm1:nm2)
    resp@amJust                = "TopLeft"
    resp@gsnPanelFigureStringsFontHeightF = subtitle_font
    resp@gsnPanelLabelBar      = True

    gsn_panel(wks,plot(nm1:nm2),(/3,2/),resp)
end do
end

