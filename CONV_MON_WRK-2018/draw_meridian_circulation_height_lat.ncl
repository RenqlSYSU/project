;***************************************************************************
;filename:      draw_zonal_circulation_height_lon.ncl
;1 read data in (ncase,12nmonth*100years,nlat,nlon),U,omega
;2 calc last 30years average and zonal average, then draw the ts plot
;                                            by Ql_Ren
;                                           2018/07/21
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
multiple = 1000
draw_var = (/"stream_func"/) + multiple 
unit     = (/"10~S~8~N~ kg/s"/);this can be used as right string 
var_name = (/"U","V","OMEGA"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)

;area: east Asia
lats = -90 ;-20
latn = 90  ;60
lonl = 0   ;20
lonr = 360 ;220
lat_spacing = 20
lon_spacing = 30

;the path of the file
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fig_name = path + "/fig/" + draw_var + "_height_lat_" + lats + "-" + latn
fig_out  = "pdf";"x11"

f    = addfile(filename(0) + var_name(0) + ".1-100.nc","r")
vars = f->$var_name(0)$(0,:,{lats:latn},{lonl:lonr})
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
nlev   = dimsizes(vars&level)
ncase  = dimsizes(case)
nvar   = dimsizes(var_name)
nyear  = 30 
nmonth = 12
nd_s = nmonth*100-nmonth*nyear
nd_e = nmonth*100-1
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
var   = new((/ncase,nvar,nyear*nmonth,nlev,nlat,nlon/),float)
var1  = new((/ncase,nvar,nlev,nlon /),float)
sf  = new((/ncase,nlev,nlat,nlon /),float)
vp  = new((/ncase,nlev,nlat,nlon /),float)
plot  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
do nv = 0,nvar-1,1
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,nv,:,:,:,:) = f->$var_name(nv)$(nd_s:nd_e,:,{lats:latn},{lonl:lonr})
end do
end do
;var(:,2,:,:,:,:) = - var(:,2,:,:,:,:)*multiple  ;magnify OMEGA 100 times
copy_VarMeta(vars,var(0,0,0,:,:,:))
printVarSummary(var)
;var1 = dim_avg_n_Wrap(dim_avg_n_Wrap(var,2),3)
var   := dim_avg_n_Wrap(var,2)
uv2sfvpf(var(:,0,:,:,:),var(:,1,:,:,:),sf,vp) ;get (ncase,nlev,nlat,nlon)
sf = sf/100000000
copy_VarMeta(vars,sf(0,:,:,:))
sf    := dim_avg_n_Wrap(sf,3)
printVarSummary(sf)

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
   mp_thick       = 2.0
   font_index     = 22
   vcref_font     = 0.01
   vc_thick       = 2.0

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18") 
    ;blue and red, has 18 colors,2-10 is blue
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
    res@vpHeightF        = 0.25
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = True 
    res@tiXAxisString   = "Longitude"
    ;res@tiYAxisString   = ""
    
    ;res@tmXBMode        = "Explicit"
    ;res@tmXBValues      = ispan(nd_s+1,nd_e+1,12)
    ;res@tmXBLabels      = ispan(ny_e-nys+1,ny_e,1)
    ;res@tmXBValues      = ispan(1,12,1)
    ;res@tmXBLabels      = months

    ;res@tmYLMode        = "Manual"
    ;res@tmYLTickStartF  = lats(na)
    ;res@tmYLTickEndF    = latn(na)
    ;res@tmYLTickSpacingF= lat_spacing

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
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    
    res@tmXMajorGrid                   = True
    res@tmXMajorGridLineDashPattern    = 2
    res@tmXMajorGridLineColor          = "black"
    res@tmXMajorGridThicknessF         = 2.0
    res@tmYMajorGrid                   = False
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "black"
    res@tmYMajorGridThicknessF         = 2.0
    
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
    
    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of vector--------------
    res@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    res@vcVectorDrawOrder          = "PostDraw"
    res@vcLineArrowThicknessF      = vc_thick
    res@vcMonoLineArrowColor       = True
    res@vcLineArrowColor           = "black"
    res@vcMinDistanceF             = 0.025
    ;res@vcRefMagnitudeF            = 10.0             ; define vector ref mag
    res@vcRefLengthF               = 0.05           ; define length of vec ref
    res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    res@vcRefAnnoOrthogonalPosF    = -0.5           ; move ref vector along Y-axis
    res@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    res@vcRefAnnoPerimOn           = False         ;close the border
    res@vcRefAnnoFont              = font_index
    res@vcRefAnnoFontHeightF       = vcref_font
    res@vcRefAnnoString1On         = False
    res@vcRefAnnoString2On         = True
    ;res@vcRefAnnoString2           = "10m/s"
    res@vcLineArrowHeadMinSizeF   = 0.01
    res@vcLineArrowHeadMaxSizeF   = 0.01

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString  = "U&OMEGA*" + multiple  
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Horizontal"
    resp@lbLabelStride    = 2
    resp@pmLabelBarWidthF  = 0.5
    resp@pmLabelBarHeightF = 0.05
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

vcref = 60
minlvl = (/-1.6/)
maxlvl = (/ 1.6/)
spacig = (/ 0.2/)

do nc = 0, ncase-1,1
    res@gsnLeftString   = pre_case(nc) + ") " + case(nc) + " (" + lats + "-" + latn + "N)"
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig        ; set contour interval
    res@cnMaxLevelValF        := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl       ; set min contour level
    res@vcRefMagnitudeF       = vcref             ; define vector ref mag
    res@vcRefAnnoString2      = vcref+"m/s"
    ;plot(nc) = gsn_csm_pres_hgt_vector(wks,var1(nc,1,:,:),var1(nc,0,:,:),var1(nc,2,:,:),res)
    ;print(max(var1(nc,1,:,:)))
    ;print(min(var1(nc,1,:,:)))
    plot(nc) = gsn_csm_pres_hgt(wks,sf(nc,:,:),res)
    print(max(sf(nc,:,:)))
    print(min(sf(nc,:,:)))
end do
    gsn_panel(wks,plot,(/ncase,1/),resp)
end

