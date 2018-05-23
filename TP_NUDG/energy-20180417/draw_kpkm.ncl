;*********************************************************************
;filename:      draw_12p3X1_month_preci-2c_shad.ncl
;input_data:    Month_ave_year_varience-6h.nc 
;Draw the contour map of monthly variance of 27 years 
;And there are three cases,1AMIP-CTRL,2TP-NUDG,3 is the difference between the above two
;the third case will do the significance testing, the dot mean the difference is significance
;The data(3variables,3cases,12month,3levs) are all save in one file
;                                                     by Ql Ren
;                                                    2018/01/08
;*********************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h"/) + "-Kp_Km_trans_10avg_vert_int1979.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/kpkm"
    case = (/"CTRL","NG6h","NG6h - CTRL" /)
fig_out  = "pdf";"x11"
unit     = "J"
ncase = dimsizes(case) 
draw_var = (/"KpKm"/)  ;this is the variable used to draw, can used as center string
var_name = (/"kpkm"/)  ;this is the variable in file 

;area: east Asia
lats = -20     
latn = 60      
lonl = 30  
lonr = 150
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case)
nvar   = dimsizes(draw_var)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var  = new((/ncase,nlat,nlon/),double)
vari = new((/ncase,nlat,nlon/),float)
prob = new((/ncase,nlat,nlon/),float)

;the number of the panel plot
    plot   = new((/ncase/),graphic)
    line   = new((/ncase/),graphic)
    gtop   = new((/ncase/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-2, 1
    f = addfile(filename(nf),"r")
    var(nf,:,:)  = f->$var_name$(0,{lats:latn},{lonl:lonr}) ;get an ctrl (ncase,nlev,nlat,nlon)
end do
var(2,:,:) = var(1,:,:) - var(0,:,:)
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"

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
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 1          ; set min contour level
    ;res@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 65                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-1

;----------Set res sig-----------
    res_sig                  = True       ; use plot options
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@gsnFrame         = False
    res_sig@gsnDraw          = False
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnLineThicknessF = 2
    ;res_sig@cnLineColor     = "black"
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 1.
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

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
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString   = draw_var

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    ;resp@txFontHeightF = title_font 
    ;resp@txFontThicknessF = title_thick

    do nc = 0, ncase-1,1
        res@gsnLeftString = case(nc)
        ;if(nc.eq.2) then
        ;    res@gsnRightString  = "siglvl: " + siglvl
            ;res@cnMinLevelValF        = -5 ;-20         ; set min contour level
            ;res@cnMaxLevelValF        = 5  ;20          ; set max contour level
            ;res@cnLevelSpacingF       = 1          ; set contour interval
        ;end if
        plot(nc) = gsn_csm_contour_map_ce(wks,var(nc,:,:),res)
        line(nc) = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
        print(max(var(nc,{lats:latn},{lonl:lonr})))
        print(min(var(nc,{lats:latn},{lonl:lonr})))
        gtop(nc) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(nc),gtop(nc))
    end do 

    gsn_panel(wks,plot,(/3,1/),resp)
end

