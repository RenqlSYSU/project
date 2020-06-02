;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
draw_var = (/"EKE","uwnd"/)  ;this is the variable used to draw, can used as center string
unit     = (/"J"  ,"m/s" /)    ;this can be used as right string 
;var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;months   = (/"DJF","MAM","JJA","SON"/)
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
case = (/"NCEP1","ERA-Interim"/)
zonal = False ; True means to calc anomaly based on zonal average
na = 0  ; "TP","EA","NH","CA"
nv = 0
nc = 1

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5 

;area: east Asia
area = (/"TP","EA","NH","CA"/)
lats = (/ 35 , 15 , 0  , 15 /);20 ;-20
latn = (/ 35 , 55 , 89 , 55 /) ;55 ;60
lonl = (/ 30 ,110 , 1  ,0   /) ;30 ;20
lonr = (/150 ,130 ,359 ,70  /) ;150;220
lev  = (/1000,925,850,700,600,500,400,300,200,150,100/)

do na = 0,0,1
;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127"
filename = path + "/mdata/ERA-Interim_Clim_month_ave_" + draw_var(nv) + "-" + cb + "-" + ca + "day.nc"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
;gtopfile = "/home/ys17-19/data-observation/NCEP_NCAR/hgt.sfc.nc"
;filename = path + "/mdata/" + case(nc) + "_Clim_month_ave_" + draw_var(nv) + ".nc" 
fig_name = path + "/fig/1p1x4_vertical_section_" + area(na) + draw_var(nv)
;fig_title  = "Climate " + draw_var(nv) + " (" + unit(0) + ")"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

nmonth = dimsizes(months)

;the number of the panel plot
    plot     = new((/nmonth/),graphic)

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
var := f->$var_name(0)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
var := dim_avg_n_Wrap(var,2)
printVarSummary(var)
if(zonal) then 
var_ave := f->$var_name(0)$(:,:,{lats(na):latn(na)},:) ;get an array (ncase,nmonth,nlat,nlon)
var = var - dim_avg_n_Wrap(var_ave,3)
end if

printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
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

    ;res@cnFillPalette = "precip3_16lev"
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
    res@tmXBTickSpacingF   = 30       ;10
    res@tmXBTickStartF     = lonl(na) ;lats(na)
    res@tmXBTickEndF       = lonr(na) ;latn(na)

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

if(nv.eq.0) then
    minlvl = (/ 3 /);
    maxlvl = (/ 48/);
    spacig = (/ 3 /);
end if

if(nv.eq.1) then
if(zonal) then 
    minlvl = (/-21/);
    maxlvl = (/ 21/);
    spacig = (/ 3 /);
else
    minlvl = (/-42/);
    maxlvl = (/ 42/);
    spacig = (/ 6 /);
end if
end if

;    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nm = 0,nmonth-1,1
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig        ; set contour interval
    res@cnMaxLevelValF        := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl       ; set min contour level
    res@gsnLeftString   = pre_case(nm) + ") " + case(nc) + " " + months(nm)
    ;res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_var(nv) + " (" + unit(nv) + ")"
    res@gsnRightString  = lats(na) + "-" + latn(na) + "E " + draw_var(nv) + " (" + unit(nv) + ")"
    plot(nm) = gsn_csm_pres_hgt(wks,var(nm,:,:),res)
    print(max(var(nm,:,:)))
    print(min(var(nm,:,:)))
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/3,4/),resp)
;    gsn_panel(wks,plot,(/1,4/),resp)
end do 
end

