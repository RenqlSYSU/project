;*********************************************************************
;filename:      draw_1p3X1_JJApreci_teleconnect.ncl
;input_data:    -Areas_JJApreci_telec.nc
;Draw the Asia teleconnection with the JJA precipitation index in four areas
;And there are there case,1AMIP-CTRL,2TP-NUDG24h,3 TP-NUDG6h
;                                            by Ql_Ren
;                                           2018/01/24
;*********************************************************************************
begin
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-TPJJAvort_daily_telec_avg.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/" + "1p1x3_TPJJAvort_daily_telec_avg"
    case = (/"CTRL","NG24h","NG6h"/)
pre_case = (/"(a","(b","(c"/)
fig_out  = "pdf";"x11"
;fig_out  = "eps";"x11"
unit     = ""
area     = (/"CTP"/) ;here used as the right title
;the data in file is corr(area,lat,lon),the area is India,Tibet,Philippines,Maritime_continent
var_name = "corr"
nsample  = 92

;plot area 
lats = -15
latn = 55
lonl = 30
lonr = 150
siglvl = 0.05

;define the dimension of the data
f = addfile(gtopfile,"r")
phis  = f->PHIS({lats:latn},{lonl:lonr})
phis  = phis/9.8   ;convert unit from m2/s2 to gpm
nlat  = dimsizes(phis&lat)
nlon  = dimsizes(phis&lon)
ncase = dimsizes(case) 
narea = dimsizes(area)
plot  = new((/ncase*narea/),graphic)
plot2 = new((/ncase*narea/),graphic)
line  = new((/ncase*narea/),graphic)
line2 = new((/ncase*narea/),graphic)
gtop  = new((/ncase*narea/),graphic)
var   = new((/ncase,narea,nlat,nlon/), float) ;the variable of teleconnection data
prob  = new((/ncase,narea,nlat,nlon/), float) ;the variable of teleconnection data

;area of TP ,PS and WMC
minlat = (/25 , 30/)
maxlat = (/40 , 35/)
minlon = (/70 , 80/)
maxlon = (/105,100/)
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   tick_thick     = 8.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

;================================================================
;read data
;====================================================================
do nf = 0, ncase-1, 1
    f             = addfile(filename(nf),"r")
    var(nf,0,:,:) = f->$var_name$(0,{lats:latn},{lonl:lonr})
end do
prob = rtest(var,nsample,0)
prob = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
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
    res@gsnMajorLatSpacing = 15
    res@gsnMinorLatSpacing = 15
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;----------------------------------------------------------------------
;set the form of the coordinate
;------------------------------------------------------------------------
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

    res@lbLabelBarOn         = False   ; the label bar
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
    
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.7
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "grey0"
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    
    res_sig@gsnCenterString  = ""
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 4.0
    
    lnres2  = lnres
    lnres2@gsLineDashPattern    = 16
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
    res_gtop@gsnCenterString  = ""
    res_gtop@gsnLeftString    = ""
    res_gtop@gsnRightString   = ""

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Horizontal";"Vertical" 
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.08
    resp@pmLabelBarOrthogonalPosF     = -0.005
    resp@pmLabelBarParallelPosF       = 0
    resp@lbLabelStride    = 1
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = tick_font-0.01
    resp@lbLabelFontThicknessF = tick_thick

do nc = 0,ncase-1,1
do na = 0,narea-1,1
    nn = na +1
    np = narea*nc + na
    res@gsnRightString   = area(na) 
    res@gsnLeftString    = pre_case(nc) + ") " + case(nc)
    plot(np)   = gsn_csm_contour_map_ce(wks,var(nc,na,:,:),res)
    plot2(np)  = gsn_csm_contour(wks,prob(nc,na,:,:),res_sig)
    overlay(plot(np),plot2(np))
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    line(np)  = gsn_add_polyline( wks, plot(np), lonx(:,0), latx(:,0), lnres )
    line2(np) = gsn_add_polyline( wks, plot(np), lonx(:,na+1), latx(:,na+1), lnres2 )
    print(max(var(nc,na,{lats:latn},{lonl:lonr})))
    print(min(var(nc,na,{lats:latn},{lonl:lonr})))
end do
end do 
    gsn_panel(wks,plot,(/narea,ncase/),resp)
end

