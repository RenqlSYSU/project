;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
year_option = True ;calc yearly DJF 
var_option = -5  ;0 means Qd; 1 means f_Qd; 2 means dzdt; 3 means duTdx; 4 means dzdt induced by b1
jet_option = True  ;20191205, control the draw of jet stream
nl = 1  ;(/925,850,700,500,200/)
na = 0 ; "NH","ME","EA","NA"

months = (/"DJF"/)
lev    = (/925,850,700,500,200/) 
case   = "ERA-Interim"
;case   = "NCEP1"

if(var_option.eq.-5) then 
    draw_var = (/"dTdy"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dTdy"/)  ;use to read variables
    unit     = (/"10~S~-5~N~K/m"/)  ;this can be used as right string 
end if
if(var_option.eq.-4) then 
    draw_var = (/"intEKE","hghEKE","lowEKE"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"intEKE","hghEKE","lowEKE"/)  ;use to read variables
    unit     = (/"10~S~4~N~J","10~S~4~N~J","10~S~4~N~J"/)  ;this can be used as right string 
end if
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
area = (/"NH","ME","EA","NA"/)
lats = (/ 10 , 10 , 10 , 10 /);10 ;-20
latn = (/ 60 , 89 , 89 , 89 /) ;55 ;60
lonl = (/ 0  , 0  ,105 , 250/) ;30 ;20
lonr = (/360 , 70 ,175 , 320/) ;150;220
density = 1.1 ;km/m3
g = 9.8 ;m2/s2

;the path of the file
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
;gtopfile = "/home/ys17-19/data-observation/NCEP_NCAR/hgt.sfc.nc"
if(var_option.le.-4)
    filename = path + "data/" + case + "_DJF_intEKE-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.le.-2.and.var_option.ge.-3)
    filename = path + "data/" + case + "_DJF_forc_te-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.le.0.and.var_option.ge.-1)
    filename = path + "data/" + case + "_DJF_forc_ex-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.eq.1)
    filename = path + "data/" + case + "_DJF_dzdt_ex-year.nc"
    fig_name = path + "fig/1p2x3_horizontal_" + var_name(nnv)
end if
if(var_option.ge.2)
    filename = path + "data/" + case + "_DJF_dzdt_te-year.nc"
    fig_name = path + "fig/1p2x3_horizontal_" + var_name(nnv)
end if

pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
fig_out    ="pdf";"x11"

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)})) ;hgt for NCEP1,unit is gpm ;z for ERA-I, unit is m2/s2
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nlat   = dimsizes(phis&latitude )  ;
nlon   = dimsizes(phis&longitude)  ;
nvar   = dimsizes(draw_var)
var_m  = new((/nvar,nlat,nlon/),double)

;the number of the panel plot
    plot     = new((/nvar/),graphic)
    line     = new((/nvar/),graphic)
    gtop     = new((/nvar/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
if(jet_option) then
    f = addfile(path + "data/" + case + "_DJF_uwnd-year.nc","r")
    u = f->season_ave(:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)})
    u_m = dim_avg_n_Wrap(u,0)
    u_m = mask(u_m, u_m.lt.20 , False) 
    printVarSummary(u_m)
end if

if(var_option.eq.-5) then
    f = addfile(path + "data/" + case + "_DJF_air-year.nc","r")
    t = dim_avg_n_Wrap(f->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    
    opt    = 0     ;used by center_finite_diff_n, no meanging 
    loncyc = True ;used by center_finite_diff_n
    cyclic = False ;used by center_finite_diff_n
    ;means NOT treated as cyclic in differenced coordinates
    ;the end values will use a one-sided difference scheme for the end points. 
    a  = 6378388 ;the radius of earth, m
    pi = atan(1.0)*4
    lat = t&lat
    lat = lat*pi/180.0
    coslat = conform(var_m(0,:,:),cos(lat),0)
    var_m(0,:,:) = tofloat(center_finite_diff_n(t*coslat,lat,cyclic,opt,0)/a/coslat)
    var_m = var_m*100000
    copy_VarMeta(t,var_m(0,:,:))
end if

if(var_option.le.0.and.var_option.ge.-4)
f = addfile(filename,"r")
if(var_option.eq.-4) then 
    do nv = 0, nvar-1,1
    if(year_option) then 
    var_m(nv,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m(nv,:,:)  = f->$var_name(nv)$({lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    end do
else
    do nv = 0, nvar-1,1
    if(year_option) then 
    var_m(nv,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m(nv,:,:)  = f->$var_name(nv)$({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    end do
end if
end if
if(var_option.ge.1)
f = addfile(filename,"r")
    if(year_option) then 
    var_m  = dim_avg_n_Wrap(f->$var_name(nnv)$(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m  = f->$var_name(nnv)$(select,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end if
printVarSummary(var_m)

if(var_option.eq.-4) then 
    var_m = var_m/10000
end if

if(var_option.le.-2.and.var_option.ge.-3)
    var_m(1,:,:) = var_m(1,:,:)*10^11
    var_m(0,:,:) = var_m(0,:,:)*24*3600
end if

if(var_option.eq.-1) then 
    var_m(2,:,:)   = var_m(2,:,:)*10^11
    var_m(0:1,:,:) = var_m(0:1,:,:)*24*3600
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

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    if(var_option.eq.-4) then 
    gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    else
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    end if
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.6
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 150.
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = 330 ;lonr(na)
    res@mpMinLonF        = -30 ;lonl(na)
    res@gsnMajorLatSpacing = 20
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
    res@pmLabelBarOrthogonalPosF     = 0.2
    res@pmLabelBarParallelPosF       = 0.5
    ;res@pmLabelBarWidthF  = 0.7
    ;res@pmLabelBarHeightF = 0.1

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    
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

;---------set the attribute of the contour-----------------
    res_line = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = True;False
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True ;False
    res_line@cnLineLabelFontHeightF = 0.005
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 5
    
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait" ;
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

if(var_option.eq.-5) then
    minlvl = (/-7 /)
    maxlvl = (/ 7 /)
    spacig = (/ 1 /)
end if
if(var_option.eq.-4) then
    minlvl = (/2  ,2  ,2  /);-8,-8,-8,-8,
    maxlvl = (/130,130,130/); 8, 8, 8, 8,
    spacig = (/8  ,8  ,8  /); 1, 1, 1, 1,
end if
if(var_option.le.-2.and.var_option.ge.-3)
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
    if(var_option.gt.-4) then 
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
    end if
do nv = 0, nvar-1,1
    np = nv
    if(var_option.gt.0) then 
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa " + draw_var(nv)
    res@gsnRightString  = var_name(nnv) + " (" + unit(nnv) + ")"
    end if
    if(var_option.le.0) then 
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa" 
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")"
    end if
    if(var_option.eq.-4) then 
    res@gsnLeftString   = pre_case(np) + ") 200hPa uwnd"
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")"
    end if
    ;res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np) = gsn_csm_contour_map_ce(wks,var_m(nv,:,:),res)
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    print(max(var_m(nv,:,:)))
    print(min(var_m(nv,:,:)))

    if(jet_option) then 
    line(np) = gsn_csm_contour(wks,u_m,res_line)
    overlay(plot(np),line(np))
    end if
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/nvar,1/),resp)
end

