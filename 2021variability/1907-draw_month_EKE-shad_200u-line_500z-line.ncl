;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
draw_var = (/"preci" ,"preci-8-80day"/) ;_filter
unit     = (/"mm/day","mm~S~2~N~/day~S~2~N~"/)     
draw_varname =(/"preci" ,"preci-8-80day"/)     
lev = (/925,500,300,200/)
nl = 1
nv = 0
jet_option  = 0 ;0 = no jet, 1 = jet axis, 2 = jet contour
vect_option = False;True ;
cont_option = False;True ;
cont_var  = (/"hgt"/)  ;this is the variable used to draw, can used as center string
cont_unit = (/"hgt"/)    ;this can be used as right string 
nvc = 0 
barshare = False;True ;
barvert  = False;True ;

gtopfile= "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
path    = "/home/ys17-19/renql/project/2021variability"
case    = (/"CN05.1","APHRO","ERAI","NCEP1","NCEP2"/) 
nc  = 0
nc2 = 2
pre_case = (/"(a","(b","(c","(d","(a","(b","(g","(h","(i","(j","(k","(l"/)
fig_out  = "pdf";"x11";"eps";

if(nv.le.1) then 
    fig_name = path + "/fig/" + case(nc) + "-month_" + draw_var(nv)
else
    fig_name = path + "/fig/" + case(nc) + "-month_" + lev(nl) + "_" + draw_var(nv)
    if(vect_option) then 
    fig_name = path + "/fig/" + case(nc) + "-month_" + lev(nl) + "wind_" + draw_var(nv)
    end if
    if(cont_option) then 
    fig_name = path + "/fig/" + case(nc) + "-month_" + lev(nl) + draw_var(nv) + "_" + cont_var(nvc)
    end if
end if

;area: east Asia
lats = (/-89,14 ,15 /);-20
latn = (/ 89,55 ,55 /);60
lonl = (/ 0 ,69 ,90 /);20
lonr = (/360,140,140/);220
siglvl = 0.05
lat_spacing = 10
lon_spacing = 20
na = 2

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

;the number of the panel plot
    plot   = new((/24/),graphic)
    plot2  = new((/24/),graphic)
    wind   = new((/24/),graphic)
    gtop   = new((/24/),graphic)
    axis   = new((/24/),graphic)
    jetl   = new((/24/),graphic)

;define the thickness and height of the title,tick and label
if(na.eq.0) then 
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   vcref_font     = 0.01
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 1.5
   font_index     = 22
else
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.0
   mp_thick       = 2.5
   font_index     = 22
end if

np = -1
wks = gsn_open_wks(fig_out,fig_name)
do nc = 0, 1, 1
var_option = 0
if(var_option.eq.0) then 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
nnv = 0
imonth = ispan(0,3,1)
filename = path + "/data/" + case(nc) + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
end if
if(var_option.eq.1) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nnv = 0
imonth = ispan(0,11,1);(/4,5,6,7/)
filename = path + "/data/" + case(nc) + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
end if

;================================================================
;read data
;====================================================================
do nv = 0, 1, 1
if(jet_option.ge.1) then 
axlon = phis&longitude
axlat = new((/nmonth,nlon/),double)
lat1 = phis&latitude
;print(lat1)
;print(ind(lat1.eq.20))
u   = new((/nmonth,nlat,nlon/),float)
u_f = addfile(filename+"uwnd.nc","r")
u   = u_f->$var_name(0)$(imonth,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon):,:);
do nm = 0, nmonth-1,1
do inln = 0, nlon-1,1
    lat_ind = maxind(u(nm,{10:48},inln))
    axlat(nm,inln) = lat1(lat_ind)
end do
end do
printVarSummary(axlat)
end if

if(vect_option) then 
f  = addfile(path+"/mdata/"+case(nc2)+"_clim_month_ave_uwnd.nc","r")
u := f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
f  = addfile(path+"/mdata/"+case(nc2)+"_clim_month_ave_vwnd.nc","r")
v := f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if

f = addfile(filename + draw_var(nv) + ".nc","r")
if(nv.le.1) then
shad  := f->$var_name(nnv)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
else
shad  := f->$var_name(nnv)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if

if(cont_option) then
f = addfile(filename+cont_var(nvc) + ".nc","r")
if(nvc.gt.1) then
cont  := f->$var_name(0)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
else
cont  := f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if
end if
printVarSummary(shad)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================

;do na = 0, 5, 1
print("draw lats = " + lats(na) + ", latn = " + latn(na) + ", lonl = " + lonl(na) + ", lonr = " + lonr(na))
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "landscape" ;"portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
if(na.eq.2.or.na.eq.5) then
    res@mpCenterLonF     = 0.
else
    res@mpCenterLonF     = 180.
end if
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
    res@gsnAddCyclic     = False
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
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 1          ; set min contour level
    ;res@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
   
;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
if(na.eq.0) then 
    resvc@vcMinDistanceF            = 0.01
    resvc@vcRefLengthF              = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
else
    resvc@vcMinDistanceF            = 0.015 ;0.025
    resvc@vcRefLengthF              = 0.05 ;0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.01;0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.01;0.015
end if
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
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
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor           = "grey"

;----------Set res line for contour-----------
    res_line                  = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False       ; Fill contours
    res_line@cnLinesOn        = True      ; turn off contour lines
    res_line@cnLineLabelsOn   = False      ; turn off contour line labels
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.0
    res_line@cnMonoLineColor  = True ;False
    if(jet_option.ge.1) then
    res_line@cnLineColor      = "dimgray"
    else
    res_line@cnLineColor      = "black"
    end if
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 2
    res_line@gsnContourPosLineDashPattern = 0
    res_line@cnLevelSelectionMode  = "AutomaticLevels"   ; set manual contour levels
    ;res_line@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res_line@cnLevelSpacingF       = 20 ;        ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;----------Set res line-----------
    res_jet                  = True
    res_jet@gsnDraw          = False
    res_jet@gsnFrame         = False
    res_jet@gsnAddCyclic     = False
    res_jet@cnFillOn         = False       ; Fill contours
    res_jet@cnLinesOn        = True      ; turn off contour lines
    res_jet@cnLineLabelsOn   = False      ; turn off contour line labels
    res_jet@cnInfoLabelOn    = False
    res_jet@cnLineThicknessF = 4
    res_jet@cnMonoLineColor  = True ;False
    res_jet@cnLineColor      = "blueviolet"
    if(jet_option.eq.1) then 
    res_jet@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_jet@cnLevels         = (/20/)
    end if
    if(jet_option.eq.2) then 
    u = mask(u, u.lt.10 , False) 
    res_jet@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_jet@cnLevelSpacingF      = 10        ; set contour interval
    end if
    res_jet@lbLabelBarOn         = False   ; the label bar
    res_jet@gsnCenterString      = ""
    res_jet@gsnLeftString    = " "
    res_jet@gsnRightString   = " "

    res_axis                      = True
    res_axis@gsLineColor          = "black"
    res_axis@gsLineDashPattern    = 1
    res_axis@gsLineThicknessF     = 3.5
    ;res_axis@gsLineLabelString    = "Westly"

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
    res_gtop@cnLevels         = (/3000/)
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

    resp  = True    
    resp@gsnMaximize   = True
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "landscaoe";"portrait"
    resp@gsnPanelLabelBar = True
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.8
    else
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.8 ;0.05
    resp@pmLabelBarHeightF = 0.05;0.8
    end if
    else
    res@lbLabelBarOn         = True   ; the label bar
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarZone               = 0
    if(barvert) then
    res@lbOrientation        = "Vertical";
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.17
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.05
    res@pmLabelBarOrthogonalPosF     = 0.2 ;y direction
    res@pmLabelBarParallelPosF       = 0.5
    end if
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(na.eq.0) then 
    vcref         = (/10, \  ;the level of 850hpa    ctrl,nudg,diff
                      10, \  ;the level of 700hpa    ctrl,nudg,diff
                      7 , \  ;the level of 500hpa    ctrl,nudg,diff
                      20/)   ;the level of 200hpa    ctrl,nudg,diff
else
    vcref         = (/7 , \  ;the level of 850hpa    ctrl,nudg,diff
                      7 , \  ;the level of 700hpa    ctrl,nudg,diff
                      7 , \  ;the level of 500hpa    ctrl,nudg,diff
                      15/)   ;the level of 200hpa    ctrl,nudg,diff
end if

if(nvc.eq.0) then ;Z3
    cont_spacing = (/20,30,50,80/)
end if

minlvl = (/1 ,0 /)  ;0.13 0.5 2  ;,-80
maxlvl = (/16,75/)  ;0.58 3.5 62 ;, 80
spacig = (/1 ,5 /)  ;0.03 0.2 4  ;, 10

if(nv.eq.2) then
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
   ; res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
   res@cnFillColors  := (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
else
    gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    res@cnFillColors  := (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)
    ;res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
end if

;do nm = 0,nmonth-1,1
do nm = 2,2,1;0,nmonth-1,1
    np = np + 1
    res@gsnLeftString   = pre_case(np) + ") " + case(nc) + " " + months(imonth(nm))
    if(cont_option) then 
    res@gsnRightString  = lev(nl) + "hPa " + draw_varname(nv) + " (" + unit(nv) + ") & " + cont_var(nvc) + " (" + cont_unit(nvc) + ") "
    else
    if(nv.le.1) then
    res@gsnRightString  = "";draw_varname(nv) + " (" + unit(nv) + ")" 
    else
    res@gsnRightString  = lev(nl) + "hPa " + draw_varname(nv) + " (" + unit(nv) + ") "
    end if
    end if
    res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np)    = gsn_csm_contour_map_ce(wks,shad(nm,:,:),res)

    if(vect_option) then 
    resvc@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl)+"m/s"
    resvc2@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc2@vcRefAnnoString2     = vcref(nl)+"m/s"
    wind(np)    = gsn_csm_vector(wks,u(nm,:,:),v(nm,:,:),resvc)
    overlay(plot(np),wind(np))
    end if
    
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    if(cont_option)
    res_line@cnLevelSpacingF    = cont_spacing(nl) 
    plot2(np)   = gsn_csm_contour(wks,cont(nm,:,:),res_line)
    overlay(plot(np),plot2(np))
    end if

    if(jet_option.ge.1) then 
    if(jet_option.eq.1) then 
    axis(np) = gsn_add_polyline( wks, plot(np), axlon, axlat(nm,:), res_axis )
    end if
    jetl(np) = gsn_csm_contour(wks,u(nm,:,:),res_jet)
    overlay(plot(np),jetl(np))
    end if
    print(max(shad(nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    print(min(shad(nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    ;draw(plot(np))
    ;frame(wks)
end do 

;gsn_panel(wks,plot,(/nmonth,1/),resp)
;gsn_panel(wks,plot,(/6,2/),resp)
;gsn_panel(wks,plot,(/5,1/),resp)
end do
end do
;resp@txString      = draw_varname(nv) + " (" + unit(nv) + ")"
gsn_panel(wks,plot,(/2,2/),resp)
end

