;*********************************************************************
;filename:      draw_JJA_qtran_div.ncl
;input_data:    ERA-Interim_Clim_daily_qtran.nc 
;                                                     by Ql Ren
;                                                    2019/03/10
;*********************************************************************************
begin
draw_var = (/"climate","steady","transient"/) 
var_name_u = (/"qu","ste_qu","edd_qu"/)
var_name_v = (/"qv","ste_qv","edd_qv"/)
unit     = "kg/(s*m)" ;"K"   ;this can be used as right string 

;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG"
path1  = "/home/ys17-19/data-observation/ERA-interim/"
filename = path + "/observation-20190127/mdata/" + "ERA-Interim_Clim_daily_qtran.nc"
gtopfile = path1 + "gtopo_1.5x1.5.nc"
fig_name = path + "/observation-20190127/fig/ERA-Interim_JJA_qtran_shad"
;fig_title= "Moisture Flux Divergence (10~S~-5~N~ kg/s/m~S~2~N~)"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i"/)
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
fig_out    ="pdf";"png";"x11";

;area: east Asia
lats = (/15 ,-15/) ;-20
latn = (/55 ,55 /) ;60
lonl = (/70 ,30 /) ;20
lonr = (/150,150/) ;220
na = 0

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)}))
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
printVarSummary(phis)
nmonth = 1 
nvar   = dimsizes(draw_var)
nlat   = dimsizes(phis&latitude)
nlon   = dimsizes(phis&longitude)
var_u  = new((/nmonth,nvar,nlat,nlon/),float)
var_v  = new((/nmonth,nvar,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/nmonth*nvar/),graphic)
    plotvc   = new((/nmonth*nvar/),graphic)
    line     = new((/nmonth*nvar/),graphic)
    gtop     = new((/nmonth*nvar/),graphic)

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
f = addfile(filename,"r")
time = cd_calendar(f->time,0)
select_time := ind((time(:,1).le.8).and.(time(:,1).ge.6))

do nm = 0,nmonth-1,1
;select_time := ind(time(:,1).eq.(nm+5))
do nv = 0,nvar-1,1
    print("read " + var_name_u(nv))
    print("read " + var_name_v(nv))
    var_u(nm,nv,:,:) = dim_avg_n_Wrap(f->$var_name_u(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}),0) 
    var_v(nm,nv,:,:) = dim_avg_n_Wrap(f->$var_name_v(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}),0) 
end do
end do
var  = sqrt(var_u^2 + var_v^2)
;var = uv2dv_cfd ( var_u, var_v, var_u&lat, var_v&lon, 2)*100000
copy_VarMeta(var_u,var)

printVarSummary(var_u)
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.5
   vc_thick       = 2.5
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"GMT_polar") 
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
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = 15
    res@gsnMinorLatSpacing = 15
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

    res@gsnCenterString      = ""
    res@gsnLeftString      = ""
    res@gsnRightString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    
;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 1          ; set min contour level
    ;res@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick
    res@pmLabelBarOrthogonalPosF = 0.2
    res@pmLabelBarParallelPosF   = 0.5

;----------set the attribute of the vector-------------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMinDistanceF             = 0.035
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = -0.78           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.25            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@vcLineArrowHeadMinSizeF   = 0.01
    resvc@vcLineArrowHeadMaxSizeF   = 0.01

    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

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
    res_gtop@cnLineThicknessF = mp_thick
    res_gtop@cnLineColor      = "black"
    res_gtop@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_gtop@cnLevels         = (/1500/)
    res_gtop@gsnCenterString      = ""
    res_gtop@gsnLeftString    = " "
    res_gtop@gsnRightString   = " "

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"  ;"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    vcref  = (/200,200,20/)  ;the level of 850hpa    ctrl,nudg,diff

if(na.eq.1) then
    maxlvl = (/24 ,24 , 8 /) ;used for East Asian moist transiation
    minlvl = (/-24,-24,-8 /)
    spacig = (/ 3 ,3  , 1 /)
else
    minlvl = (/30 ,30 , 0/) ;used for East Asian moist transiation
    maxlvl = (/180,180,30/)
    spacig = (/10 ,10 , 2/)
end if
   
    ;res@cnFillPalette = "BlueDarkRed18"
    res@cnFillPalette = "precip3_16lev"
    do nm = 0, nmonth-1,1
    do nv = 0 ,nvar-1,1
        np = 3*nm + nv
        res@gsnLeftString   = pre_case(nv) + ") " + "JJA" 
        res@gsnRightString  = draw_var(nv) ;months(nm+4)
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF      := spacig(nv)        ; set contour interval
        res@cnMaxLevelValF       := maxlvl(nv) ;32    ; set max contour level
        res@cnMinLevelValF       := minlvl(nv)       ; set min contour level
        resvc@vcRefMagnitudeF   = vcref(nv)                ; define vector ref mag
        resvc@vcRefAnnoString2  = vcref(nv) + " " + unit
        plot(np)   = gsn_csm_contour_map_ce(wks,var(nm,nv,:,:),res) 
        plotvc(np) = gsn_csm_vector(wks,var_u(nm,nv,:,:),var_v(nm,nv,:,:),resvc)
        gtop(np)   = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),plotvc(np))
        overlay(plot(np),gtop(np))
        ;line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var_u(nm,nv,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        print(min(var_u(nm,nv,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        print(max(var(nm,nv,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        print(min(var(nm,nv,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    end do
    end do
;resp@txString      = fig_title
gsn_panel(wks,plot,(/nmonth,nvar/),resp)
end

