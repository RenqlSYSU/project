;*********************************************************************
;    Draw the regression patterns for 850 hPa wind regressed by the JJA precipitation index in four areas
;use the gray shaded to sign the regression coefficient of U is through the significance test
;and the spot means the coefficient of V is through the test
;    There are there case,1AMIP-CTRL,2TP-NUDG,3 is the difference between the above two
;The data(3cases,4areas) are all save in one file,and are all named "var" in the script.
;                                                               by QlRen
;                                                               2017/12/02
;*********************************************************************************

begin
;the path of the file
filename   = "/users/yangsong3/Model_Team/F/clm-month/Areas_preci_JJA_reg850wind.nc"
fig_name   = "/users/yangsong3/Model_Team/F/clm-month/fig/4panel_3X1_JJA_preci_reg850wind"
fig_title  = "JJA preci index regression patterns"
left_title = (/"AMIP_CTRL","TP_NUDG","AMIP_CTRL - TP_NUDG"/)
fig_out    = "pdf";"x11"
unit       = "m/s"
area      = (/"India","Tibet","Philippines","Maritime_continent"/) ;here used as the right title
ncase = 3 
narea = dimsizes(area)
var_name = (/"corr","tval"/)
;plot area 
lats = -15
latn = 55
lonl = 30
lonr = 150
siglvl = 0.05
;define the dimension of the data
    nlat  = 85
    nlon  = 113
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0

;the number of the panel plot
    plot   = new((/3/),graphic)
    var    = new((/ncase,narea,2,nlat,nlon/), float) ;the variable of teleconnection data,2 means U and V
    tlat   = new((/ncase,narea,2/), float) ;store the lat which through the significance test
    tlon   = new((/ncase,narea,2/), float) ;store the lon which through the significance test

;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
var  = f->$var_name(0)$
tval = f->$var_name(1)$
print(var(1,1,1,20,50))

;==============================================================
;draw figure in twelve months (one month is a fgure)
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

;--------------------------------------------------------
;set the attribute of the map
;-------------------------------------------------------------
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
    res@cnMinLevelValF        = -0.9          ; set min contour level
    res@cnMaxLevelValF        = 0.9        ; set max contour level
    res@cnLevelSpacingF       = 0.1        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
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
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font

do na = 0,narea-1,1
    res@gsnRightString       = area(na) 
    do nc = 0, ncase-1,1
        res@gsnLeftString    = left_title(nc)
        plot(nc) = gsn_csm_contour_map_ce(wks,var(nc,na,:,:),res)
        print(max(var(nc,na,{lats:latn},{lonl:lonr})))
        print(min(var(nc,na,{lats:latn},{lonl:lonr})))
    end do
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title 
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    gsn_panel(wks,plot,(/3,1/),resp)
end do 
end

