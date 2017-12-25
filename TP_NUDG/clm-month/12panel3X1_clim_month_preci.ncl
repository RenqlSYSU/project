;**********************************************************************
;draw the monthly climate precipitation in the  AMIP-CTRL and TP-NUDG
;**********************************************************************

begin
;===========================================================
;variable list
;=========================================================
;define the dimension of the data
    ntime = 12
    nlev  = 1
    nlat  = 192
    nlon  = 288
    nvar  = 3 
;area: east Asia
    lats = -20
    latn = 60
    lonl = 20
    lonr = 220
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0

;the filename of the read-in-data and figure, as well as the title of the figure
    filein     = "/users/yangsong3/Model_Team/F/clm-month/" + "preci.nc"
    fig_name   = "/users/yangsong3/Model_Team/F/clm-month/fig/" +\
                 "east_asia_12panel3X1_clim_month_preci_diff"
    fig_title  = "Climate Month Preci"
    left_title = (/"AMIP_CTRL","TP_NUDG","AMIP_CTRL - TP_NUDG"/)
    fig_out    = "pdf";"x11"_
    unit       = "mm/day"
;the number of the panel plot
    plot   = new((/3/),graphic)
    var    = new((/nvar,ntime,nlat,nlon/), float)
    var_name   = (/"CTRLpreci","NUDGpreci","DIFFpreci"/)
    months = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nev","Dec"/)

;================================================================
;read data
;====================================================================
    f = addfile(filein,"r")
do i = 0,nvar-1,1  
    var(i,:,:,:) = f->$var_name(i)$  
end do
print(var(1,11,20,50))

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
wks = gsn_open_wks(fig_out,fig_name)
do m = 0,11,1
; m=1
   ; wks = gsn_open_wks(fig_out,fig_name+months(m))
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True
    res@gsnPaperOrientation = "portrait"

;--------------------------------------------------------
;set the attribute of the map
;-------------------------------------------------------------
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

;----------------------------------------------------------------------
;set the form of the coordinate
;------------------------------------------------------------------------
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

;--------------------------------------------------------------------------
;set the attribute of the contour
;--------------------------------------------------------------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 4          ; set min contour level
    res@cnMaxLevelValF        = 36          ; set max contour level
    res@cnLevelSpacingF       = 4          ; set contour interval

    res@gsnSpreadColors      = True               ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnRightString       = unit
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    
    print(m)
    do i = 0, nvar-1,1
        res@gsnLeftString    = left_title(i)
        if(i.eq.2) then
            res@gsnSpreadColors = False
            res@cnMinLevelValF        = -5         ; set min contour level
            res@cnMaxLevelValF        = 5          ; set max contour level
            res@cnLevelSpacingF       = 1          ; set contour interval
        end if
        plot(i) = gsn_csm_contour_map_ce(wks,var(i,m,:,:),res)
    print(max(var(i,m,{lats:latn},{lonl:lonr})))
    print(min(var(i,m,{lats:latn},{lonl:lonr})))
    end do
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title + " in " + months(m)
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    gsn_panel(wks,plot,(/3,1/),resp)
;draw(plot)
;frame(wks)
end do 
end

