;*********************************************************************************
;Draw the first three EOF spatial patterns and principal components of region uwnd 
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;                                            by Ql_Ren
;                                           2020/02/23
;*********************************************************************************
begin
area = (/"ME","EA","NA"/);,"GB"
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 0  ,105 , 250/);, 0   ;30 ;20
lonr = (/ 70 ,175 , 320/);, 360 ;150;220
lon_spacing = 20
lat_spacing = 10
narea = dimsizes(area)
neval = 4

path   = "/home/ys17-19/renql/project/2020MiddleEastJet"
filename   = path + "/data/eof_u_6kinds-" + area + ".nc"
fig_name   = path + "/fig/eof_u_pattern-" + area 
fig_out    = "pdf" ;"eps";"x11"
pre_case   = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)
plot   = new((/6/),graphic)
line   = new((/6/),graphic)

;==============================================================
;draw figure 
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   font_index     = 22

do na = 0, narea-1,1
    wks = gsn_open_wks(fig_out,fig_name(na))
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True

;------------set the form of the coordinate-----------------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick

    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font

    resxy = res
    resxy@vpWidthF  = 0.4
    resxy@vpHeightF = 0.3
    resxy@gsnStringFontHeightF = subtitle_font/2
    resxy@tmXBLabelFontHeightF = tick_font/2
    resxy@tmYLLabelFontHeightF = tick_font/2
    
    resxy@xyDashPatterns    = (/0,0/) ;,0,0
    resxy@xyLineThicknesses = (/4.0, 4.0/) ;, 4.0, 4.0  ; make second line thicker
    resxy@xyLineColors      = (/"black","red"/); ,"darkorange3","chartreuse4"
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = -0.5 ;floor(min(var)); -0.4;
    ;resxy@trYMaxF =  0.5 ;ceil(max(var)) ; 0.2 ;
    resxy@tiXAxisOn = True
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    resxy@gsnYRefLine           = 0.0
    resxy@gsnYRefLineColor      = "black"
    resxy@gsnYRefLineThicknessF = 5.0
    resxy@gsnYRefLineDashPattern= 0

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
;    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -0.08      ; set min contour level
    res@cnMaxLevelValF        = 0.08       ; set max contour level
    res@cnLevelSpacingF       = 0.01        ; set contour interval

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    ;res@pmLabelBarZone               = 0
    ;res@pmLabelBarOrthogonalPosF     = 0.2
    ;res@pmLabelBarParallelPosF       = 0.5

    res_vert = res
    res_vert@vpHeightF        = 0.6
    res_vert@vpWidthF         = 0.8
    res_vert@gsnPresHgtHeightLabelOn = False
    res_vert@tiYAxisOn = False
    res_vert@tmYRLabelsOn = False
    res_vert@tmXBMode = "Manual"

;------------set the attribute of the map----------------
    res@gsnAddCyclic     = False
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;---------set the attribute of the contour-----------------
    res_line = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True ;False
    res_line@cnLineLabelFontHeightF = 0.015
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    res_line@gsnContourZeroLineThicknessF = 8
    levels = new((/10/),integer)
    levels(0) = 0
    levels(1:9) = ispan(5,85,10)
    res_line@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    res_line@cnLevels             = levels       ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    resp@gsnPaperOrientation = "landscape"

f2 = addfile(path + "/data/ERA-Interim_DJF_uwnd-year.nc","r")
vars := dim_avg_n_Wrap(f2->season_ave(:,:,{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
f3 = addfile(path + "/data/ERA-Interim_DJF_ps-year.nc","r")
ps   := dim_avg_n_Wrap(f3->season_ave(:,{lats(na):latn(na)},{lonl(na):lonr(na)}),0)

f      = addfile(filename(na),"r")
nyear  = dimsizes(f->eof_ts(0,0,:))
eof0 := f->eof0
eof1 := f->eof1
eof2 := f->eof2
eof3 := f->eof3
eof4 := f->eof4
eof5 := f->eof5

;----------horizontal-------------------------
nk = 0
np = 0
sig  = eofunc_north(eof0@pcvar,nyear,True)
subtitle_color = where(sig,"red","black") 
do ne = 0,neval-2,1
    res@gsnLeftString    = area(na) + nk 
    res@gsnRightStringFontColor   = subtitle_color(ne)
    res@gsnRightString            = eof0@eval(ne)+"  "+eof0@pcvar(ne)+"%" 
    plot(np) = gsn_csm_contour_map_ce(wks,eof0(ne,:,:),res)
    line(np) = gsn_csm_contour(wks,vars({200},:,:),res_line)
    overlay(plot(np),line(np))
    np = np + 1
end do

nk = 1
ptop = 70 
iopt = 0 ;Set to zero. Currently not used.
dp  := dpres_plevel(vars&lev*100,ps,ptop,iopt)
var1:= dim_sum_n(vars*dp,0)/dim_sum_n(dp,0);vertical average
copy_VarMeta(vars(0,:,:),var1)
printVarSummary(var1)
sig  = eofunc_north(eof1@pcvar,nyear,True)
subtitle_color = where(sig,"red","black") 
do ne = 0,neval-2,1
    res@gsnLeftString    = area(na) + nk 
    res@gsnRightStringFontColor   = subtitle_color(ne)
    res@gsnRightString            = eof1@eval(ne)+"  "+eof1@pcvar(ne)+"%" 
    plot(np) = gsn_csm_contour_map_ce(wks,eof1(ne,:,:),res)
    line(np) = gsn_csm_contour(wks,var1,res_line)
    overlay(plot(np),line(np))
    np = np + 1
end do
gsn_panel(wks,plot,(/2,3/),resp)

;----------vertical section-----------------------
    res_vert@tmXBTickSpacingF   = lat_spacing 
    res_vert@tmXBTickStartF     = lats(na)
    res_vert@tmXBTickEndF       = latn(na)
    res_line@gsnPresHgtHeightLabelOn = False
    res_line@tiYAxisOn = False

nk = 2
np = 0
var2 := dim_avg_n_Wrap(vars,2)
printVarSummary(var2)
sig  = eofunc_north(eof2@pcvar,nyear,True)
subtitle_color = where(sig,"red","black") 
do ne = 0,neval-2,1
    res_vert@gsnLeftString    = area(na) + nk 
    res_vert@gsnRightStringFontColor   = subtitle_color(ne)
    res_vert@gsnRightString            = eof2@eval(ne)+"  "+eof2@pcvar(ne)+"%" 
    plot(np) = gsn_csm_pres_hgt(wks,eof2(ne,:,:),res_vert)
    line(np) = gsn_csm_pres_hgt(wks,var2,res_line)
    overlay(plot(np),line(np))
    np = np + 1
end do

    res_vert@tmXBTickSpacingF   = lon_spacing 
    res_vert@tmXBTickStartF     = lonl(na)
    res_vert@tmXBTickEndF       = lonr(na)
nk = 4
var4 := dim_avg_n_Wrap(vars,1)
printVarSummary(var4)
sig  = eofunc_north(eof4@pcvar,nyear,True)
subtitle_color = where(sig,"red","black") 
do ne = 0,neval-2,1
    res_vert@gsnLeftString    = area(na) + nk 
    res_vert@gsnRightStringFontColor   = subtitle_color(ne)
    res_vert@gsnRightString            = eof4@eval(ne)+"  "+eof4@pcvar(ne)+"%" 
    plot(np) = gsn_csm_pres_hgt(wks,eof4(ne,:,:),res_vert)
    line(np) = gsn_csm_pres_hgt(wks,var4,res_line)
    overlay(plot(np),line(np))
    np = np + 1
end do
gsn_panel(wks,plot,(/2,3/),resp)

;----------xy plot-------------------------------------
    resxy@tiXAxisString = "lat"
    resxy@tmXBMode      = "Explicit"
    resxy@tmXBValues    = ispan(lats(na),latn(na),lat_spacing)
    resxy@tmXBLabels    = ispan(lats(na),latn(na),lat_spacing)
nk = 3
np = 0
var3 = dim_avg_n_Wrap(var1,1)
var3 = var3/200
printVarSummary(var3)
sig  = eofunc_north(eof3@pcvar,nyear,True)
subtitle_color = where(sig,"red","black") 
do ne = 0,neval-2,1
    resxy@gsnLeftString    = area(na) + nk 
    resxy@gsnRightStringFontColor   = subtitle_color(ne)
    resxy@gsnRightString            = eof3@eval(ne)+"  "+eof3@pcvar(ne)+"%" 
    plot(np) = gsn_csm_xy(wks,var3&lat,(/eof3(ne,:),var3/),resxy)
    np = np + 1
end do

    resxy@tiXAxisString = "lon"
    resxy@tmXBValues   := ispan(lonl(na),lonr(na),lon_spacing)
    resxy@tmXBLabels   := ispan(lonl(na),lonr(na),lon_spacing)
nk = 5
var5 = dim_avg_n_Wrap(var1,0)
var5 = var5/200
printVarSummary(var5)
sig  = eofunc_north(eof5@pcvar,nyear,True)
subtitle_color = where(sig,"red","black") 
do ne = 0,neval-2,1
    resxy@gsnLeftString    = area(na) + nk 
    resxy@gsnRightStringFontColor   = subtitle_color(ne)
    resxy@gsnRightString            = eof5@eval(ne)+"  "+eof5@pcvar(ne)+"%" 
    plot(np) = gsn_csm_xy(wks,var5&lon,(/eof5(ne,:),var5/),resxy)
    np = np + 1
end do
gsn_panel(wks,plot,(/2,3/),resp)
delete([/res,resxy,res_vert,res_line/])
end do
end
