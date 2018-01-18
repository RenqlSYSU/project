;**********************************************************************
;draw the monthly climate 500hpa height in the  AMIP-CTRL and TP-NUDG
;**********************************************************************

begin
;===========================================================
;variable list
;=========================================================
;define the dimension of the data
    ntime = 12
    lev   = (/850,500,200/)  ;the hpa of the level
    ;nlev  = dimsizes(lev)
    nlat  = 192
    nlon  = 288
    ;nvar  = 3  ;u,v,omega
    ;ncase = 3  ;CTRL,NUDG,CTRL-NUDG
;area: east Asia
    lats = -20
    latn = 60
    lonl = 20
    lonr = 220
;define the contour level
    maxlevel      = (/(/0.27,0.27,0.045/), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/0.27,0.27,0.045/), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/0.09,0.09,0.045/)/)   ;the level of 200hpa    ctrl,nudg,diff

    minlevel      = -(/(/0.27,0.27,0.045/), \  ;the level of 850hpa    ctrl,nudg,diff 
                       (/0.27,0.27,0.045/), \  ;the level of 500hpa    ctrl,nudg,diff
                       (/0.09,0.09,0.045/)/)   ;the level of 200hpa    ctrl,nudg,diff

    levelspacing  = (/(/0.03,0.03,0.005/), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/0.03,0.03,0.005/), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/0.01,0.01,0.005/)/)   ;the level of 200hpa    ctrl,nudg,diff

    vcref         = (/(/15,15,5/), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/30,30,5/), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/50,50,5/)/)   ;the level of 200hpa    ctrl,nudg,diff

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

;the filename of the read-in-data and figure, as well as the title of the figure
    path       = "/users/yangsong3/Model_Team/F/clm-month/"
    pre_name   =(/"AMIP_C5PM.cam.h1.ESM.clim.month.",\
                 "AMIP_C5PM_TP_NUDG.cam.h1.ESM.clim.month."/)
    var_name   = "Z3"
    
    fig_name   = path + "fig/" + "east_asia_12panel_3X1_clim_month_height_diff"
    fig_title  = "Climate Month Height"
    left_title = (/"AMIP_CTRL","TP_NUDG","AMIP_CTRL - TP_NUDG"/)
    right_title = "m"
    fig_out    = "pdf";"x11"

;the number of the panel plot
    plot   = new((/3/),graphic)
    var    = new((/ncase,ntime,nlat,nlon/), float)
    months = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nev","Dec"/)

h = 1
;do h = 0,dimsizes(lev)-1,1   ; the loop of the height
;================================================================
;read data
;====================================================================
   ; f = addfile(filein,"r")
do i = 0,1,1
    filein = path + pre_name(i) + var_name + ".nc"
    f = addfile(filein,"r")
    var(i,:,:,:) = f->$var_name$(:,{lev(h)},:,:)
end do
var(2,:,:,:) = var(0,:,:,:) - var(1,:,:,:)    ;the difference of the ctrl and nudg
print(var(1,11,20,50))

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
wks = gsn_open_wks(fig_out,fig_name+lev(h))
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
    res@gsnScalarContour = True
    res@cnFillOn         = False 
    res@cnLinesOn        = True
    res@cnLineLabelsOn   = True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = minlevel(h)          ; set min contour level
    ;res@cnMaxLevelValF        = maxlevel(h)          ; set max contour level
    ;res@cnLevelSpacingF       = levelspacing(h)          ; set contour interval

    res@gsnSpreadColors      = False ;True               ; turn of the full color map
    ;res@gsnSpreadColorStart  = 11                ; start at color 17
    ;res@gsnSpreadColorEnd    = 2                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1

;---------------------------------------------------------------------
;set the attribute of the vector
;---------------------------------------------------------------------------
     res@vcGlyphStyle               = "CurlyVector"
     res@vcVectorDrawOrder          = "PostDraw"
     res@vcLineArrowThicknessF      = vc_thick
     res@vcMinDistanceF             = 0.015
     ;res@vcRefMagnitudeF            = 30.0             ; define vector ref mag
     res@vcRefLengthF               = 0.05           ; define length of vec ref
     res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
     ;res@vcRefAnnoSide              = "Bottom"          
     ;res@vcRefAnnoJust              = "TopRight"
     res@vcRefAnnoOrthogonalPosF    = 0.5           ; move ref vector along Y-axis
     res@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
     res@vcRefAnnoPerimOn           = False         ;close the border
     res@vcRefAnnoFontHeightF       = vcref_font
     res@vcRefAnnoString1On         = False
     res@vcRefAnnoString2On         = True
     ;res@vcRefAnnoString2           = "30m/s"

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnRightString       = right_title
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    
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

