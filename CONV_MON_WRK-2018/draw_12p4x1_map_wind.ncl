;***************************************************************************
;filename:      draw_regin_map_preci.ncl
;1 read data in (ncase,12nmonth*100years,nlat,nlon),PRECC + PRECL
;2 calc last 30years monthly average, then draw the ts plot
;                                            by Ql_Ren
;                                           2018/07/19
;******************************************************************************
begin
draw_var = (/"wind"/)
unit     = (/"m/s"/)    ;this can be used as right string 
var_name = (/"U","V"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)

;area: east Asia
lats = -90 ;-20
latn = 90  ;60
lonl = 0  ;20
lonr = 360 ;220
lev  = (/850,500,200/)
nl = 0
vcref  = 10
lat_spacing = 30
lon_spacing = 40

;the path of the file
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fig_name = path + "/fig/animation/" + draw_var + "_map_" + lats + "-" + latn
fig_out  = "png";"pdf";"x11"

f    = addfile(filename(0) + var_name(0) + ".1-100.nc","r")
vars = f->$var_name(0)$(0,{lev(nl)},{lats:latn},{lonl:lonr})
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
ncase  = dimsizes(case)
nvar   = dimsizes(var_name)
nyear  = 30 
nmonth = 12
nd_s = nmonth*100-nmonth*nyear
nd_e = nmonth*100-1
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
var   = new((/ncase,nvar,nyear*nmonth,nlat,nlon/),float)
var1  = new((/ncase,nvar,nmonth      ,nlat,nlon/),float)
plot  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
do nv = 0,nvar-1,1
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,nv,:,:,:) =  f->$var_name(nv)$(nd_s:nd_e,{lev(nl)},{lats:latn},{lonl:lonr}) ;+ means downward to earth
end do
end do
copy_VarMeta(vars,var(0,0,0,:,:))
printVarSummary(var)

do nm = 0,11,1
var1(:,:,nm,:,:) = dim_avg_n_Wrap(var(:,:,nm:(nyear*nmonth-1):12,:,:),2)
end do

;do ns = 0,nyear-1,1
;    nms = ns*nmonth
;    nme = ns*nmonth + 11 
;    var1(:,ns,:,:) = dim_avg_n(var(:,nms:nme,:,:),1) ;calc the annual mean
;end do
;var1 := dim_avg_n(var1(:,69:99,:,:),1)
;copy_VarMeta(vars,var1(0,:,:))
;printVarSummary(var1)

;===============================================
;draw the figure
;=============================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   vcref_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   vc_thick       = 2
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    ;gsn_define_colormap(wks,"BlueDarkRed18") 
    ;;blue and red, has 18 colors,2-10 is blue
    gsn_define_colormap(wks,"precip3_16lev") 
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

;--------set the attribute of vector--------------
    res@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    res@vcVectorDrawOrder          = "PostDraw"
    res@vcLineArrowThicknessF      = vc_thick
    res@vcMonoLineArrowColor       = True
    res@vcLineArrowColor           = "black"
    res@vcMinDistanceF             = 0.03
    res@vcRefMagnitudeF            = vcref             ; define vector ref mag
    res@vcRefAnnoString2           = vcref + " m/s"
    res@vcRefLengthF               = 0.05           ; define length of vec ref
    res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;res@vcRefAnnoSide              = "Bottom"          
    ;res@vcRefAnnoJust              = "TopRight"
    res@vcRefAnnoOrthogonalPosF    = 0.4           ; move ref vector along Y-axis
    res@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    res@vcRefAnnoPerimOn           = False         ;close the border
    res@vcRefAnnoFont              = font_index
    res@vcRefAnnoFontHeightF       = vcref_font
    res@vcRefAnnoString1On         = False
    res@vcRefAnnoString2On         = True
    res@vcLineArrowHeadMinSizeF   = 0.012
    res@vcLineArrowHeadMaxSizeF   = 0.012

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@gsnPanelLabelBar = False
    resp@pmLabelBarWidthF  = 0.1
    resp@pmLabelBarHeightF = 0.8
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@lbOrientation = "Vertical"
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do nm = 0,12-1,1 
;ny = floor((nt)/12)+1   ;used draw every month every year
;print(ny)
;nm = toint(nt-(ny-1)*12)
;print(nm)
do nc = 0, ncase-1,1
    nv = 0
    res@gsnLeftString   = pre_case(nc) + ")" + case(nc)
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ") " + months(nm) ;lev(nl) + "hPa " +  
    plot(nc)   = gsn_csm_vector_map(wks,var1(nc,0,nm,:,:),var1(nc,1,nm,:,:),res)
    print(max(var(nc,0,nm,:,:)))
    print(min(var(nc,0,nm,:,:)))
    ;plot(nc)   = gsn_csm_contour_map(wks,var1(nc,:,:),res)
    ;print(max(var1(nc,:,:)))
    ;print(min(var1(nc,:,:)))
end do
;    resp@txString      = "Year" + ny
    gsn_panel(wks,plot,(/2,2/),resp)
end do
end

