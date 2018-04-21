;*********************************************************************
;filename:      draw_4p3X1_JJApreci_teleconnect.ncl
;input_data:    -Areas_JJApreci_telec.nc
;Draw the Asia teleconnection with the JJA precipitation index in four areas
;And there are there case,1AMIP-CTRL,2TP-NUDG24h,3 TP-NUDG6h
;                                            by Ql_Ren
;                                           2018/01/24
;*********************************************************************************
begin
path       = "/users/yangsong3/renql/project/TP_NUDG"
;path       = "/home/ys17-19/renql/project/TP_NUDG" 
filename   = path + "/analysis/mdata/" + (/"CTRL","NG6h"/) + "-Areas_JJApreci_daily_telec.nc"
fig_name   = path + "/phy_panel-20180328/fig/" + "p2x2_JJApreci_telec"
    case   = (/"CTRL","NG6h"/)
fig_out    = "pdf";"x11"
unit       = ""
area       = (/"India","Tibet","Philippines","Maritime_continent"/) ;here used as the right title
var_name   = "corr"
nsample    = 27*92

;plot area 
lats = -15
latn = 55
lonl = 30
lonr = 150
siglvl = 0.05

;define the dimension of the data
f     = addfile(filename(0),"r")
vars  = f->$var_name$(0,{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
ncase = dimsizes(case) 
narea = dimsizes(area)
plot  = new((/ncase/),graphic)
plot2 = new((/ncase/),graphic)
line  = new((/ncase/),graphic)
line2 = new((/ncase/),graphic)
var   = new((/ncase,narea,nlat,nlon/), float) ;the variable of teleconnection data
prob  = new((/ncase,narea,nlat,nlon/), float) ;the variable of teleconnection data

minlat = (/25 ,20 ,25  ,5   ,-10 /)
maxlat = (/40 ,35 ,35  ,20  ,7.5 /)
minlon = (/70 ,60 ,90  ,120 ,85  /)
maxlon = (/105,80 ,103 ,140 ,120 /)
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0

;================================================================
;read data
;====================================================================
do nf = 0, ncase-1, 1
    f             = addfile(filename(nf),"r")
    var(nf,:,:,:) = f->$var_name$(:,{lats:latn},{lonl:lonr})
end do
prob = rtest(var,nsample,0)
prob = mask(prob,prob.lt.siglvl,True)
printVarSummary(var)
printVarSummary(prob)

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

;----------Set res sig-----------
    res_sig                  = True
    res_sig@gsnDraw          = False
    res_sig@gsnFrame         = False
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnInfoLabelOn    = False
    ;res_sig@cnLineThicknessF = 2
    ;res_sig@cnLineColor     = "black"
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 1.
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/0.01/)
    ;res_sig@cnFillPattern        = "GMT_gray"
    ;res_sig@cnFillColors         = (/"grey30","grey50"/) ;
    
    res_sig@gsnCenterString      = ""
    res_sig@gsnRightStringFontHeightF = subtitle_font
    res_sig@gsnLeftStringFontHeightF  = subtitle_font
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 4.0
    ;lnres@gsLineLabelString    = "TP"

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font

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

do na = 0,narea-1,1
    res@gsnRightString       = area(na) 
    do nc = 0, ncase-1,1
        res@gsnLeftString    = case(nc)
        plot(nc)   = gsn_csm_contour_map_ce(wks,var(nc,na,:,:),res)
        plot2(nc)  = gsn_csm_contour(wks,prob(nc,na,:,:),res_sig)
        overlay(plot(nc),plot2(nc))
        line(nc)  = gsn_add_polyline( wks, plot(nc), lonx(:,0), latx(:,0), lnres )
        line2(nc) = gsn_add_polyline( wks, plot(nc), lonx(:,na+1), latx(:,na+1), lnres )
        print(max(var(nc,na,{lats:latn},{lonl:lonr})))
        print(min(var(nc,na,{lats:latn},{lonl:lonr})))
    end do
    
    resp@txString      = fig_title 
    gsn_panel(wks,plot,(/3,1/),resp)
end do 
end

