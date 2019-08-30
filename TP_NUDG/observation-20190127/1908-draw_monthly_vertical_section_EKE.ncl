;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
draw_var = (/"EKE"/)  ;this is the variable used to draw, can used as center string
unit     = (/"J"/)    ;this can be used as right string 
;var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;months   = (/"DJF","MAM","JJA","SON"/)
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nv = 0

;parameter for bw_bandpass_filter
ca = 7.0 ;day
cb = 2.5 

;area: east Asia
lats = -15 ;20 ;-20
latn = 55  ;55 ;60
lonl =30  ;30 ;20  210;
lonr =150 ;150;220 330;
lat = 40 ;used to plot vertical section

;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127"
path1  = "/home/ys17-19/data-observation/ERA-interim/"
;filename = path + "/mdata/ERA-Interim_Clim_month_ave_" + draw_var(nv) + "-" + cb + "-" + ca + "day.nc"
filename = path + "/mdata/ERA-Interim_Clim_month_ave_" + draw_var(nv) + ".nc" 
gtopfile = path1 + "gtopo_1.5x1.5.nc"
fig_name = path + "/fig/1p1x4_ERA-Interim_vertical_section_" + lat + draw_var(nv)
fig_title  = "Climate " + draw_var(nv) + " (" + unit(0) + ")"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

nmonth = dimsizes(months)

;the number of the panel plot
    plot     = new((/nmonth/),graphic)
    line1    = new((/nmonth/),graphic)
    line2    = new((/nmonth/),graphic)
    gtop     = new((/nmonth/),graphic)
    axis     = new((/nmonth/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.5
   font_index     = 22

;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
var = f->$var_name(0)$(:,:,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
;print(var(1,8,{lat},{30:40}))
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
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
    res@vpHeightF        = 0.6
    res@vpWidthF         = 0.8

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
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
;    res@cnLineThicknessF = 1.0

    res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;res@cnFillPalette = "GMT_gray" ;has 10 colors 
    res@cnLevelSelectionMode  = "ManualLevels" ;"AutomaticLevels" ; set manual contour levels
    res@cnMinLevelValF        = 0    ; set min contour level
    res@cnMaxLevelValF        = 160   ; set max contour level
    res@cnLevelSpacingF       = 10  ; set contour interval
    
    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0
    res@pmLabelBarParallelPosF       = 0.5

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    
    res@gsnPresHgtHeightLabelOn = False
    res@tiYAxisOn = False
    res@tmYRLabelsOn = False
    res@tmXBMode = "Manual"
    res@tmXBTickSpacingF   = 30
    res@tmXBTickStartF     = lonl
    res@tmXBTickEndF       = lonr

;---------set the attribute of the contour-----------------
    res_line = res
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.0
    res_line@cnLineColor      = "black"
   
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do nm = 0,nmonth-1,1
    res@gsnLeftString   = pre_case(nm) + ") Obs " + months(nm)
    res@gsnRightString  = lat + "N " + draw_var(nv) + " (" + unit(nv) + ")"
    res_line@gsnLeftString   = pre_case(nm) + ") Obs " + months(nm)
    res_line@gsnRightString  = lat + "N " + draw_var(nv) + " (" + unit(nv) + ")"
    plot(nm) = gsn_csm_pres_hgt(wks,var(nm,:,{lat},:),res)
    ;plot(nm) = gsn_csm_pres_hgt(wks,var(nm,:,{lat},:),res_line)
    print(max(var(nm,:,{lat},{lonl:lonr})))
    print(min(var(nm,:,{lat},{lonl:lonr})))
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/3,4/),resp)
end

