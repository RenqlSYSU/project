;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
year_option = True ;calc yearly DJF 
var_option = 1  ;0 means Qd; 1 means f_Qd; 2 means dzdt; 3 means duTdx; 4 means dzdt induced by b1
jet_option = True  ;20191205, control the draw of jet stream
stream_option = True
months = (/"DJF"/)
lev    = (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,175,150,125, 100/)
case   = "ERA-Interim"
;case   = "NCEP1"

if(var_option.eq.-3) then 
    draw_var = (/"Qeddh","Ah"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qeddh","Ah"/)  ;use to read variables
    unit     = (/"K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.-2) then 
    draw_var = (/"Qeddl","Al"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qeddl","Al"/)  ;use to read variables
    unit     = (/"K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.-1) then 
    draw_var = (/"Qd"   ,"Qeddy","A"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qd"   ,"Qeddy","A"/)  ;use to read variables
    unit     = (/"K/Day","K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.0) then 
    draw_var = (/"Qd"   ,"Qeddy","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qd"   ,"Qeddy"/)  ;use to read variables
    unit     = (/"K/Day","K/Day","K/day"/)  ;this can be used as right string 
end if
if(var_option.eq.1) then 
    draw_var = (/"Qd"  ,"Qeddy","A"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
    nnv = 1 ;used to read dzdt induced by different forcing
    select = (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.2) then 
    draw_var = (/"Qeddl","Al"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
    nnv = 1 ;used to read dzdt induced by different forcing
    select = (/1,3/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.3) then 
    draw_var = (/"Qeddh","Ah"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
    nnv = 1 ;used to read dzdt induced by different forcing
    select = (/0,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if

;area: east Asia
area = (/"ME","EA","NA"/)
lats = (/ 10 , 10 , 10 /);10 ;-20
latn = (/ 60 , 60 , 60 /) ;55 ;60
lonl = (/ 0  ,105 , 250/) ;30 ;20
lonr = (/ 70 ,175 , 320/) ;150;220
na = 0
density = 1.1 ;km/m3
g = 9.8 ;m2/s2

;the path of the file
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"
if(var_option.le.-2)
    filename = path + "data/" + case + "_DJF_Qeddy_A_freq-year.nc"
    fig_name = path + "fig/1p3x3_height_lat_" + var_name(0)
end if
if(var_option.le.0.and.var_option.ge.-1)
    filename = path + "data/" + case + "_DJF_Qd_Qeddy_A-year.nc"
    fig_name = path + "fig/1p3x3_height_lat_" + var_name(0)
end if
if(var_option.eq.1)
    filename = path + "data/" + case + "_DJF_dzdt-year.nc"
    fig_name = path + "fig/1p2x3_height_lat_" + var_name(nnv)
end if
if(var_option.ge.2)
    filename = path + "data/" + case + "_DJF_dzdt_freq-year.nc"
    fig_name = path + "fig/1p2x3_horizontal_" + var_name(nnv)
end if

pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
fig_out    ="pdf";"x11"

f = addfile(filename,"r")
if(var_option.le.0)
    vars := f->$var_name(0)$(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
else
    vars := f->$var_name(nnv)$(0,0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if
nlat   = dimsizes(vars&lat)
nlev   = dimsizes(lev)
nvar   = dimsizes(draw_var)
narea  = dimsizes(area)
var_m  = new((/narea,nvar,nlev,nlat/),double)

;the number of the panel plot
    plot     = new((/nvar*narea/),graphic)
    line     = new((/nvar*narea/),graphic)
    msfp     = new((/nvar*narea/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.035
   subtitle_font  = 0.035
   title_font     = 0.03
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
ulev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)
unlev = dimsizes(ulev)

if(jet_option) then
u_m   = new((/narea,unlev,nlat/),float)
f = addfile(path + "data/" + case + "_DJF_uwnd-year.nc","r")
do na = 0,narea-1,1
    u := f->season_ave(:,:,{lats(na):latn(na)},{lonl(na):lonr(na)})
    u := dim_avg_n_Wrap(u,3)
    u_m(na,:,:) = dim_avg_n_Wrap(u,0)
end do
u_m = mask(u_m, u_m.lt.10 , False) 
printVarSummary(u_m)
;copy_VarMeta(vars(:,:,0),u(0,:,:))
end if

if(stream_option) then
msf_m   = new((/narea,unlev,nlat/),float)
f  = addfile(path + "data/" + case + "_DJF_vwnd-year.nc","r")
f2 = addfile(path + "data/" + case + "_DJF_ps-year.nc","r")
do na = 0,narea-1,1
    v   := dim_avg_n_Wrap(f->season_ave(:,:,{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    ps  := dim_avg_n_Wrap(f2->season_ave(:, {lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    msf := zonal_mpsi_Wrap(v(::-1,:,:),v&lat,ulev(::-1)*100,ps)
    msf_m(na,:,:) = msf(::-1,:) 
end do
msf_m = msf_m/10^11
printVarSummary(msf_m)
print(max(msf_m))
print(min(msf_m))
end if

f = addfile(filename,"r")
if(var_option.lt.0)
    do nv = 0, nvar-1,1
    do na = 0, narea-1,1
    if(year_option) then 
    var  := dim_avg_n_Wrap(f->$var_name(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var  := f->$var_name(nv)$({lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    var_m(na,nv,:,:) = dim_avg_n_Wrap(var,2)
    end do
    end do
else
    do na = 0, narea-1,1
    if(year_option) then 
    var  := dim_avg_n_Wrap(f->$var_name(nnv)$(select,:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var  := f->$var_name(nnv)$(select,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    var_m(na,:,:,:) = dim_avg_n_Wrap(var,3)
    end do
end if
copy_VarMeta(vars(:,:,0),var_m(0,0,:,:))
printVarSummary(var_m)

if(var_option.lt.-1) then 
    var_m(:,1,:,:) = var_m(:,1,:,:)*10^11
    var_m(:,0,:,:) = var_m(:,0,:,:)*24*3600
end if

if(var_option.eq.-1) then 
    var_m(:,2,:,:) = var_m(:,2,:,:)*10^11
    var_m(:,0:1,:,:) = var_m(:,0:1,:,:)*24*3600
end if

if(var_option.eq.0) then 
    var_m = var_m*24*3600 ;convert units of Q from K/s to K/Day 
end if

if(var_option.ge.1) then
if(nnv.eq.0) then
    var_m = var_m*10^3
end if
if(nnv.eq.1) then
    var_m = var_m*10^6
end if
if(nnv.eq.2) then
    var_m = var_m*10^6
end if
if(nnv.eq.3) then
    var_m = var_m*10^12
end if
end if
printVarSummary(var_m)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
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

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
   
    if(var_option.ge.1) then
    res@lbLabelBarOn         = False   ; the label bar
    else
    res@lbLabelBarOn         = True   ; the label bar
    end if
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0
    res@pmLabelBarParallelPosF       = 0.5
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.1

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
    res@tmXBTickSpacingF   = 10
    res@tmXBTickStartF     = lats(0)
    res@tmXBTickEndF       = latn(0)

;---------set the attribute of the contour-----------------
    res_msf = True
    res_msf@gsnDraw          = False
    res_msf@gsnFrame         = False
    res_msf@cnFillOn         = False
    res_msf@cnLinesOn        = True
    res_msf@cnLineLabelsOn   = True ;False
    res_msf@cnInfoLabelOn    = False
    res_msf@cnLineThicknessF = 3.0
    res_msf@cnLineColor      = "blueviolet";"gray40"
    res_msf@gsnContourZeroLineThicknessF = 0
    res_msf@gsnContourNegLineDashPattern = 2
    res_msf@gsnContourPosLineDashPattern = 0
    
    res_msf@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_msf@cnLevelSpacingF       = 1.0        ; set contour interval
    res_msf@gsnPresHgtHeightLabelOn = False
    res_msf@tiYAxisOn = False
    res_msf@lbLabelBarOn     = False   ; the label bar
    res_msf@gsnCenterString  = ""
    res_msf@gsnLeftString    = " "
    res_msf@gsnRightString   = " "

;---------set the attribute of the contour-----------------
    res_line = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True ;False
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "black"
    ;res_line@gsnContourZeroLineThicknessF = 0
    ;res_line@gsnContourNegLineDashPattern = 5
    
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@gsnPresHgtHeightLabelOn = False
    res_line@tiYAxisOn = False
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    if(var_option.ge.1) then
    resp@gsnPanelLabelBar = True
    else
    resp@gsnPanelLabelBar = False
    end if
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.06
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(var_option.lt.-1) then
    minlvl = (/-2.1,-14/);-8,-8,-8,-8,
    maxlvl = (/ 2.1, 14/); 8, 8, 8, 8,
    spacig = (/ 0.3, 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.-1) then
    minlvl = (/-2.1,-2.1,-14/);-8,-8,-8,-8,
    maxlvl = (/ 2.1, 2.1, 14/); 8, 8, 8, 8,
    spacig = (/ 0.3, 0.3, 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.0) then 
    minlvl = (/-2.1,-2.1,-2.1/);
    maxlvl = (/ 2.1, 2.1, 2.1/);
    spacig = (/ 0.3, 0.3, 0.3/);
end if

if(var_option.ge.1) then
if(nnv.eq.0) then   ;dzdt
    minlvl = (/-3.5,-3.5,-3.5,-3.5/);
    maxlvl = (/ 3.5, 3.5, 3.5, 3.5/);
    spacig = (/ 0.5, 0.5, 0.5, 0.5/);
end if
if(nnv.eq.1) then ;dudt
    minlvl = (/-35,-35,-35,-35/);-8,-8,-8,-8,
    maxlvl = (/ 35, 35, 35, 35/); 8, 8, 8, 8,
    spacig = (/ 5 , 5 , 5 , 5 /); 1, 1, 1, 1,
end if
if(nnv.eq.2) then ;dTdt
    minlvl = (/-16,-16,-16,-16/);
    maxlvl = (/ 16, 16, 16, 16/);
    spacig = (/  2,  2,  2,  2/);
end if
if(nnv.eq.3) then ;dbdt
    minlvl = (/-24,-24,-24,-24/);
    maxlvl = (/ 24, 24, 24, 24/);
    spacig = (/  3,  3,  3,  3/);
end if
end if

    ;res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do na = 0, narea-1,1
do nv = 0, nvar-1,1
    np = na*nvar + nv
    res@gsnLeftString   = pre_case(np) + ") " + area(na) 
    if(var_option.ge.1) then 
    res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_var(nv) + " (" + unit(nnv) + ")"
    else
    res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_var(nv) + " (" + unit(nv) + ")"
    end if
    ;res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np) = gsn_csm_pres_hgt(wks,var_m(na,nv,:,:),res)
    print(max(var_m(nv,na,:,:)))
    print(min(var_m(nv,na,:,:)))

    if(jet_option) then 
    line(np) = gsn_csm_pres_hgt(wks,u_m(na,:,:),res_line)
    overlay(plot(np),line(np))
    end if
    
    if(stream_option) then 
    msfp(np) = gsn_csm_pres_hgt(wks,msf_m(na,:,:),res_msf)
    overlay(plot(np),msfp(np))
    end if
end do 
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/narea,nvar/),resp)
end

