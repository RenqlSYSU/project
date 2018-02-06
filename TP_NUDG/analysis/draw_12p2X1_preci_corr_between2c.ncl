;*********************************************************************
;filename:      draw_12p2X1_preci_corr_between2c.ncl
;input_data:    CTRL-Preci_corr_with_ctrl.nc 
;Draw the correlation shaded map between the CTRL precipitation and NUDG 
;to see whether the precipitation is in cocurrent phase.
;                                                     by Ql Ren
;                                                    2018/01/20
;*********************************************************************************
begin
;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
;path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/analysis/mdata/" + (/"NUDG24h","NUDG6h"/) + "-Preci_corr_with_ctrl.nc"
fig_name   = path + "/analysis/fig/" + (/"12p2X1_preci_corr_between2c_month","4p2X1_preci_corr_between2c_season"/)
fig_title  = "Preci Correlation with CTRL in " ;+season
    case   = (/"NUDG24h","NUDG6h"/)
fig_out    ="pdf";"x11"
unit       = "mm/day"   ;this can be used as right string 
draw_var = (/"preci"/)  ;this is the variable used to draw, can used as center string
var_name = (/"month_corr","season_corr"/) ;this is the var_name in data file

mon_season = 0 ;0 is to draw monthly data,1 is to draw seasonal data
if (mon_season.eq.0) then
    months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;    nsample  = (/31,28,31,30,31,30,31,31,30,31,30,31/)*27
else
    months   = (/"DJF","MAM","JJA","SON"/)
;    nsample  = (/90,92,92,91/)*27
end if
nsample = 27
print(nsample)

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
lev  = (/850,500,200/)  ;can used as right string
siglvl = 0.05

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
nmonth = dimsizes(months)
nlat   = getfilevardimsizes(f,"lat")
nlon   = getfilevardimsizes(f,"lon")
var  = new((/ncase  ,nmonth,nlat,nlon/),float)
vari = new((/ncase  ,nmonth,nlat,nlon/),float)
prob = new((/ncase  ,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase/),graphic)
    plot2    = new((/ncase/),graphic)
    line     = new((/ncase/),graphic)

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
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf),"r")
    var(nf,:,:,:)  = f->$var_name(mon_season)$ ;get an array (ncase,nmonth,nlat,nlon)
end do

do nf = 0 ,ncase-1,1
do nm = 0 ,nmonth-1,1
prob(nf,nm,:,:) = rtest(var(nf,nm,:,:),nsample,0)
;prob(nf,nm,:,:) = rtest(var(nf,nm,:,:),nsample(nm),0)
end do
end do

prob = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
printVarSummary(var)
printVarSummary(prob)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name(mon_season))
    gsn_define_colormap(wks,"GMT_polar") 
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
    res@gsnMajorLonSpacing = 20
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
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

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    ;res@cnLineThicknessF = 0
    ;res@gsnContourZeroLineThicknessF = 4.0
    ;res@cnLineColor      = "blue"
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -0.27          ; set min contour level
    res@cnMaxLevelValF        = 0.27        ; set max contour level
    res@cnLevelSpacingF       = 0.03        ; set contour interval

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 2

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
    res@gsnRightString  = "siglvl: " + siglvl

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

do nm = 0,nmonth-1,1
    do nc = 0, ncase-1,1
        res@gsnLeftString   = case(nc)
        plot(nc)     = gsn_csm_contour_map_ce(wks,var(nc,nm,:,:),res)
        plot2(nc)    = gsn_csm_contour(wks,prob(nc,nm,:,:),res_sig)
        overlay(plot(nc),plot2(nc))
        line(nc)     = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
        print(max(var(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nm,{lats:latn},{lonl:lonr})))
    end do 

    resp@txString      = fig_title + months(nm) 
    gsn_panel(wks,plot,(/2,1/),resp)
end do ;for the month 
end

