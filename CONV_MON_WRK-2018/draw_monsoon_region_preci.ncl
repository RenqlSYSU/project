;***************************************************************************
;filename:      draw_monsoon_region_preci.ncl
;1 read data in (ncase,nmonth*nyears,nlat,nlon),PRECC + PRECL
;2 calc last 30years monthly average, 
;3 calc summer(5,6,7,8,9) and winter(11,12,1,2,3) and annual total average preci (SH NH seperate)
;4 calc monsoon index (summer-winter>2mm/day and summer*5/(annual*12)>0.55)
;5 draw the index
;                                            by Ql_Ren
;                                           2018/07/19
;******************************************************************************
begin
draw_var = (/"monsoon_region_preci"/)
unit     = (/"mm/day"/)    ;this can be used as right string 
var_name = (/"PRECC","PRECL"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)

;area: east Asia
lats = -60 ;-20
latn = 60  ;60
lonl = 0  ;20
lonr = 360 ;220
lat_spacing = 20
lon_spacing = 30

;the path of the file
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fig_name = path + "/fig/" + draw_var + lats + "-" + latn
fig_out  = "pdf";"x11"

f    = addfile(filename(0) + var_name(0) + ".1-100.nc","r")
vars = f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
ncase  = dimsizes(case)
nyear  = 30 
nmonth = 12
nd_s = nmonth*100-nmonth*nyear
nd_e = nmonth*100-1
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
var   = new((/ncase,nyear*nmonth,nlat,nlon/),float)
var1  = new((/ncase,nmonth      ,nlat,nlon/),float)
summer  = new((/ncase,nlat,nlon/),float)
winter  = new((/ncase,nlat,nlon/),float)
annual  = new((/ncase,nlat,nlon/),float)
index   = new((/ncase,nlat,nlon/),float)
plot  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
    nv = 0
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    fsnt =  f->$var_name(nv)$(nd_s:nd_e,{lats:latn},{lonl:lonr}) ;+ means downward to earth
    
    nv = 1
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    flnt =  f->$var_name(nv)$(nd_s:nd_e,{lats:latn},{lonl:lonr}) ;+ means upward to space

    var(nc,:,:,:) = fsnt + flnt ;+ means downward to earth
end do
var  = var*3600*24*1000 ;convert unit from m/s to mm/day
copy_VarMeta(vars,var(0,0,:,:))
copy_VarMeta(vars,summer(0,:,:))
copy_VarMeta(vars,winter(0,:,:))
printVarSummary(var)

do nm = 0,11,1
var1(:,nm,:,:) = dim_avg_n_Wrap(var(:,nm:(nyear*nmonth-1):12,:,:),1)
end do
summer(:,{lats:0},:) = dim_avg_n_Wrap(var1(:,(/10,11,0,1,2,3/),{lats:0},:),1)
summer(:,{0:latn},:) = dim_avg_n_Wrap(var1(:,4:9,{0:latn},:),1)
winter(:,{lats:0},:) = dim_avg_n_Wrap(var1(:,4:9,{lats:0},:),1)
winter(:,{0:latn},:) = dim_avg_n_Wrap(var1(:,(/10,11,0,1,2,3/),{0:latn},:),1)
annual = dim_avg_n_Wrap(var1,1)
index = mask(annual,((summer-winter).gt.2).and.((summer*6/(annual*12)).gt.0.55),True)
copy_VarMeta(vars,index(0,:,:))
printVarSummary(index)

;===============================================
;draw the figure
;=============================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    ;gsn_define_colormap(wks,"BlueDarkRed18") 
    ;;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
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
   ; res@mpProjection      = "Mollweide"       ; choose projection
   ; res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
   ; res@mpPerimOn         = False             ; turn off box around plot
   ; res@mpGridLatSpacingF = 30.               ; spacing for lat lines
   ; res@mpGridLonSpacingF = 30.               ; spacing for lon lines
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = lat_spacing 
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = True
    res@mpLandFillColor = "gray"
    res@mpFillDrawOrder = "PreDraw"
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
    
;--------set the attribute of the mark----------------
    res@cnFillOn         = True       ; Fill contours
    res@cnLinesOn        = False;True      ; turn off contour lines
    res@cnLineLabelsOn   = False      ; turn off contour line labels
    res@cnInfoLabelOn    = False
    ;res@cnLineThicknessF = 2
    ;res@cnLineColor     = "black"
    res@cnMonoFillPattern= True;False
    res@cnFillPattern    = 11
    res@cnMonoFillScale  = True
    res@cnFillScaleF     = 0.5
    res@cnMonoFillColor  = True
    res@cnFillColor      = "red"
    res@cnFillDotSizeF   = 0.5
    res@lbLabelBarOn     = False;True   ; turn off the label bar

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape"
    resp@gsnPanelLabelBar = False 
    resp@pmLabelBarWidthF  = 0.1
    resp@pmLabelBarHeightF = 0.8
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@lbOrientation = "Vertical"
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do nc = 0, ncase-1,1
    nv = 0
    res@gsnLeftString   = pre_case(nc) + ")" + case(nc)
    res@gsnRightString  = draw_var(nv)
    plot(nc)   = gsn_csm_contour_map(wks,index(nc,:,:),res)
end do
    gsn_panel(wks,plot,(/ncase,1/),resp)
end

