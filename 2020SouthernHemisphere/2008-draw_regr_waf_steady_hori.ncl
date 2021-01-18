;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
jet_option  = False;True ;;20191205, control the draw of jet stream
bar_share = False;
barvert   = False;
lev    = (/925,850,500,300,200,  100,50,10/) 
nl = 4 ;(/925,850,700,500,200/)
na = 0 ; "NH","ME","EA","NA"
DJF_option = True ;False;
siglvl= 0.05

if(DJF_option) then
months  = (/"DJF_"/)
months2 = (/"DJF"/)
else
months  = (/""/)
;months2 = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months2 = (/"DJF","MAM","JJA","SON"/)
end if
ntime = dimsizes(months2) 
case = (/"NCEP2","ERA-Interim"/)
nc = 0
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
fig_out    ="pdf";"x11"

;area: east Asia
;area = (/"NH","ME","EA","NA","TP"/)
lats = (/ 5  , 0  ,  0 ,  0 , 10 /);10 ;-20
latn = (/ 89 , 89 , 89 , 89 , 50 /) ;55 ;60
lonl = (/ 0  , 0  ,105 , 250, 40 /) ;30 ;20
lonr = (/360 , 70 ,175 , 320,130 /) ;150;220
lat_spacing = 15
lon_spacing = 20
density = 1.1 ;km/m3
g = 9.8 ;m2/s2

do var_option = 0,0,1 
if(var_option.eq.0) then 
    draw_var = (/"steady"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"wafx"  ,"wafy"  ,"edd_sf" /)  ;
    unit     = (/"m~S~2~N~/s~S~2~N~" ,"m2/s2" ,"10~S~6~N~m~S~2~N~/s"/)  ;this can be used as right string 
    shad_option = 1 ;0=divg; 1=edd_sf; 2=wafz 
end if
if(var_option.eq.1) then 
    draw_var = (/"RWS"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"rws" /)  ;
    unit     = (/"10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
end if

;the path of the file
path   = "/home/ys17-19/renql/project/2020SouthernHemisphere/"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
if(var_option.eq.0)
    filename = path + "data/" + case(nc) + "_" + months + "steady_TN01-comp.nc"
    fig_name = path + "fig/1p3x3_pcomp_horizontal_TN01_" + var_name(2)
    ;filename = path + "data/" + case(nc) + "_" + months + "steady_TN01-reg.nc"
;    fig_name = path + "fig/1p3x3_regr_horizontal_TN01_" + var_name(2)
end if
if(var_option.eq.1)
    filename = path + "data/" + case + "_" + months + "rws-reg.nc"
    fig_name = path + "fig/1p3x3_regr_horizontal_" + var_name(0)
end if

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)})) ;hgt for NCEP1,unit is gpm ;z for ERA-I, unit is m2/s2
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

f = addfile(filename,"r")
vars = f->$var_name(0)$({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
nlat   = dimsizes(vars&lat)  ;
nlon   = dimsizes(vars&lon)  ;
nvar   = dimsizes(draw_var)
if(var_option.eq.0) then 
var_m  := new((/3,ntime,nlat,nlon/),double)
prob   := new((/ntime,nlat,nlon/),double)
else
var_m  := new((/nvar,ntime,nlat,nlon/),float)
end if

;the number of the panel plot
    plot     = new((/12/),graphic)
    line     = new((/12/),graphic)
    sigl     = new((/12/),graphic)
    line1    = new((/12/),graphic)
    line2    = new((/12/),graphic)
    line3    = new((/12/),graphic)
    gtop     = new((/12/),graphic)
    vect     = new((/12/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   vcref_font     = 0.005
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue

do nl = 4, 4, 1
;================================================================
;read data
;====================================================================
if(jet_option) then
    f = addfile(path + "data/" + case + "_" + months + "uwnd-year.nc","r")
    u_m = new((/ntime,nlat,nlon/),float)
    if(DJF_option) then
    u_m(0,:,:) = dim_avg_n_Wrap(f->season_ave(:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    else
    u_m = dim_avg_n_Wrap(f->month_ave(:,:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    end if
    u_m = mask(u_m, u_m.lt.10 , False) 
    printVarSummary(u_m)
end if

f = addfile(filename,"r")
do nv = 0,2,1
    if(DJF_option) then
    var_m(nv,0,:,:) = f->$var_name(nv)$({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m(nv,:,:,:)   = f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end do

if(var_option.eq.1) then 
    var_m = var_m*10^11
    minlvl := (/-7,-3.5,-3.5/); -1.4,-1.4,-1.4-3.5,-3.5,-3.5-7,-7,-7
    maxlvl := (/ 7, 3.5, 3.5/);  1.4, 1.4, 1.4 3.5, 3.5, 3.5 7, 7, 7
    spacig := (/ 1, 0.5, 0.5/);  0.2, 0.2, 0.2 0.5, 0.5, 0.5 1, 1, 1
end if

if(var_option.eq.0) then 
;var_m = var_m*1000
vcref := (/4,150,150/) ;(/925,850,700,500,200/)
;vcref := (/150,150,150/) ;(/925,850,700,500,200/)
;vcref := (/50,50,50/) ;(/925,850,700,500,200/)
prob(0,:,:) = f->prob({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
prob   = mask(prob,prob.lt.siglvl,True)
lon_m  = conform(prob,var_m&lon,2)
prob   = mask(prob,lon_m.eq.0,False)
delete(lon_m)
if(shad_option.eq.0) then ;divg
    var_m(2,:,:,:) = uv2dv_cfd( var_m(0,:,:,:), var_m(1,:,:,:),var_m&lat, var_m&lon, 3 )*100000
    minlvl := (/-3.5,-3.5,-3.5/); -1.4,-1.4,-1.4-3.5,-3.5,-3.5-7,-7,-7
    maxlvl := (/ 3.5, 3.5, 3.5/);  1.4, 1.4, 1.4 3.5, 3.5, 3.5 7, 7, 7
    spacig := (/ 0.5, 0.5, 0.5/);  0.2, 0.2, 0.2 0.5, 0.5, 0.5 1, 1, 1
end if
if(shad_option.eq.1) then ;edd_sf
    var_m(2,0,:,:) = (1.0/1000000.0)*f->$var_name(2)$({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    minlvl := (/-7,-21,-21/);-3.5,-3.5,-3.54.2,-4.2,-4.2
    maxlvl := (/ 7, 21, 21/); 3.5, 3.5, 3.54.2, 4.2, 4.2
    spacig := (/ 1, 3 , 3 /); 0.5, 0.5, 0.50.6, 0.6, 0.6
end if
if(shad_option.eq.2) then ;wafz
    var_name(2) = "wafz"
    unit(2)     = "10~S~-8~N~m~S~2~N~/s~S~2~N~"
    var_m(2,0,:,:) = 100000000.0*f->$var_name(2)$({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    minlvl := (/-14/);-3.5,-3.5,-3.54.2,-4.2,-4.2
    maxlvl := (/ 14/); 3.5, 3.5, 3.54.2, 4.2, 4.2
    spacig := (/ 2 /); 0.5, 0.5, 0.50.6, 0.6, 0.6
end if
end if
printVarSummary(var_m)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
;    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.6
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 150.
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    if(na.eq.0) then
    res@gsnAddCyclic     = True;False
    res@mpMaxLonF        = 330 ;lonr(na)
    res@mpMinLonF        = -30 ;lonl(na)
    else
    res@gsnAddCyclic     = False;True;
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    end if
    res@gsnMajorLatSpacing = lat_spacing 
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
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
    
    ;res@lbLabelBarOn         = True   ; the label bar
    ;res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

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

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
    resvc@vcMinDistanceF             = 0.012
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
    ;resvc@vcRefMagnitudeF            = 40.0             ; define vector ref mag
    ;resvc@vcRefAnnoString2           = "40m/s"
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.3           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

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
    if(na.eq.0) then 
    res_line@gsnAddCyclic     = True;False
    else
    res_line@gsnAddCyclic     = False
    end if

;----------Set res sig-----------
    res_sig                  = True       ; use plot options
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@gsnFrame         = False
    res_sig@gsnDraw          = False
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnLineThicknessF = 2
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "gray50";"darkslategray"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "
    if(na.eq.0) then 
    res_sig@gsnAddCyclic     = True ;False;
    else
    res_sig@gsnAddCyclic     = False
    end if

;--------set the attribute of polygon--------------
    lnres1                      = True
    lnres1@gsLineColor          = "black"
    lnres1@gsLineDashPattern    = 0
    lnres1@gsLineThicknessF     = 5.0
    lnres2       = lnres1
    lnres2@gsLineDashPattern    = 2
    lnres3       = lnres1
    lnres3@gsLineDashPattern    = 16

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait" ;
    if(bar_share) then
    res@lbLabelBarOn      = False;True   ; the label bar
    resp@gsnPanelLabelBar = True
    if(barvert) then 
        resp@lbOrientation     = "Vertical"
        resp@pmLabelBarWidthF  = 0.05
        resp@pmLabelBarHeightF = 0.6
    else
        resp@lbOrientation     = "Horizontal" ;"Vertical"
        resp@pmLabelBarWidthF  = 0.8 ;0.05
        resp@pmLabelBarHeightF = 0.05;0.8
    end if
    else
    res@lbLabelBarOn      = True   ; the label bar
    resp@gsnPanelLabelBar = False
    if(barvert) then
        res@lbOrientation        = "Vertical";
        res@pmLabelBarWidthF  = 0.03
        res@pmLabelBarHeightF = 0.17
        res@pmLabelBarOrthogonalPosF     = 0.01
        res@pmLabelBarParallelPosF       = 0.5
    else
        res@lbOrientation        = "Horizontal";"Vertical";
        res@pmLabelBarWidthF  = 0.7
        res@pmLabelBarHeightF = 0.1
        res@pmLabelBarOrthogonalPosF     = 0.25 ;y direction
        res@pmLabelBarParallelPosF       = 0.5
    end if
    end if
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick


    ;res@cnFillPalette = "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
    
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)

do nt = 0, ntime-1,1
do nv = 0, nvar-1, 1
    np = nl-4
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa positive5";negitive6"
    if(var_option.eq.0) then 
    if(shad_option.eq.0) then 
    res@gsnRightString  = draw_var(nv) + " WAF  divg (" + unit(2) + ")"
    else
    res@gsnRightString  = draw_var(nv) + " WAF  " + var_name(2) + " (" + unit(2) + ")"
    end if
    else
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")"
    end if
    
    ;res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
   
    if(var_option.eq.0) then 
    plot(np) = gsn_csm_contour_map_ce(wks,var_m(2,nt,:,:),res)
    sigl(np) = gsn_csm_contour(wks,prob(nt,:,:),res_sig)
    overlay(plot(np),sigl(np))
    resvc@vcRefMagnitudeF      = vcref(nv)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nv)+unit(0)
    vect(np) = gsn_csm_vector(wks,var_m(0,nt,:,:),var_m(1,nt,:,:),resvc)
    overlay(plot(np),vect(np))
    printMinMax(var_m(0,nt,:,:),0)
    printMinMax(var_m(1,nt,:,:),0)
    printMinMax(var_m(2,nt,:,:),0)
    else
    plot(np) = gsn_csm_contour_map_ce(wks,var_m(nv,nt,:,:),res)
    printMinMax(var_m(nv,nt,:,:),0)
    end if

    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    if(jet_option) then 
    line(np) = gsn_csm_contour(wks,u_m(nt,:,:),res_line)
    overlay(plot(np),line(np))
    end if
end do 
end do 

;    resp@txString      = fig_title
;    gsn_panel(wks,plot,(/1,1/),resp)
;     gsn_panel(wks,plot,(/nvar,1/),resp)
end do 
    gsn_panel(wks,plot,(/1,1/),resp)
end do 
end

