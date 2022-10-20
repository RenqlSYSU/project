;***************************************************************************
;filename:      draw_preci_lat_time.ncl
;1 read data in (ncase,12nmonth*100years,nlat,nlon),PRECC + PRECL
;2 calc zonal mean,then plot
;                                            by Ql_Ren
;                                           2018/07/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
begin
draw_var = (/"TS"/)
unit     = (/"~S~o~N~C"/)    ;this can be used as right string 
var_name = (/"TS"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)
nv = 0
nys  = 30
nd_s = 1200-nys*12
nd_e = nd_s+nys*12-1
ny_e = (nd_e+1)/12

;area: east Asia
lats = -30 ;-20
latn = 30  ;60
lonl = 160  ;20
lonr = 260 ;220
lat_spacing = 15
lon_spacing = 30

;the path of the file
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fig_name = path + "/fig/" + draw_var(nv) + "_lat_time_" + lonl + "-" + lonr
fig_out  = "pdf";"x11"

f    = addfile(filename(0) + var_name(0) + ".1-100.nc","r")
vars = f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
ncase  = dimsizes(case)
nmonth = 12
;months = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
var   = new((/ncase,nys*12,nlat,nlon/),float)
var1  = new((/ncase,12    ,nlat,nlon/),float)
plot  = new((/ncase/),graphic)
line  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================

;---------read data and convert unit-----------------
var = 0
do nc = 0,ncase-1,1
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,:,:,:) =  f->$var_name(nv)$(nd_s:nd_e,{lats:latn},{lonl:lonr}) ;+ means downward to earth
end do
var = var-273.15 ;convert unit from k to c
printVarSummary(var)

;-------in order to draw two cycle years--------------
var12    = new((/ncase,24    ,nlat/),float)
months2  = new((/24/),string)
months2(0:11)  = months
months2(12:23) = months

;------calculate climate annual cycle-------------
do nm = 0,11,1
    var1(:,nm,:,:) = dim_avg_n_Wrap(var(:,nm:(nys*12-1):12,:,:),1)
end do
var1 := dim_avg_n_Wrap(var1,3) ;get(ncase,ntime,nlat)
;var1  = dim_rmvmean_n_Wrap(var1,1) ;romve annual mean, get(ncase,ntime,nlat)
var12(:,0:11,:)  = var1
var12(:,12:23,:) = var1
var12!0    = "case"
var12!1    = "time"
var12&time = ispan(1,24,1)
printVarSummary(var12)
var2 = var12(case|:,lat|:,time|:)

pai = 3.14
sun_point = 23.433*sin(pai/6.0*var12&time-pai/2.0)
printVarSummary(sun_point)

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
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 3.0
   mp_thick       = 3.0
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    ;gsn_define_colormap(wks,"BlueDarkRed18") 
    ;blue and red, has 18 colors,2-10 is blue
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
    res@vpHeightF        = 0.25
    res@vpWidthF         = 0.9

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiXAxisString   = "Year"
    res@tiYAxisString   = "Latitude(N)"
    
    res@tmXBMode        = "Explicit"
    ;res@tmXBValues      = ispan(nd_s+1,nd_e+1,12)
    ;res@tmXBLabels      = ispan(ny_e-nys+1,ny_e,1)
    res@tmXBValues      = ispan(1,24,1)
    res@tmXBLabels      = months2

    ;res@tmYLMode        = "Manual"
    ;res@tmYLTickStartF  = lats(na)
    ;res@tmYLTickEndF    = latn(na)
    ;res@tmYLTickSpacingF= lat_spacing

    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont         = font_index
    res@tmYLLabelFont         = font_index
    res@tmXBLabelFontHeightF     = tick_font 
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    res@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    res@tmXTOn               = False ; Turn off x-axis (top) major ticks
    res@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    res@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    res@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    
    res@tmXMajorGrid                   = True
    res@tmXMajorGridLineDashPattern    = 2
    res@tmXMajorGridLineColor          = "black"
    res@tmXMajorGridThicknessF         = 2.0
    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "black"
    res@tmYMajorGridThicknessF         = 2.0
    
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    res@cnlevelselectionmode = "automaticlevels"   ; set manual contour levels
    ;res@cnlevelselectionmode = "explicitlevels"   ; set manual contour levels
    ;res@cnlevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = -4.5          ; set min contour level
    ;res@cnMaxLevelValF        =  4.5 ;32        ; set max contour level
    ;res@cnLevelSpacingF       =  0.5        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 2
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-2

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;-----------set the attribute of the polyline---------------
    lres                = True
    lres@gsLineThicknessF = 7.0            ; line thickness   
    lres@gsLineColor      = "black"         ; color of lines
    
;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font

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

    minlvl = (/23/)
    maxlvl = (/30/)
    spacig = (/0.5/)

do nc = 0, ncase-1,1
    res@gsnLeftString   = pre_case(nc) + ")" + case(nc)
    res@gsnRightString  = lonl + "-" + lonr + "E " + draw_var + " (" + unit + ")"
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(nc)   = gsn_csm_lat_time(wks,var2(nc,:,:),res)
    line(nc)   = gsn_add_polyline(wks,plot(nc),var12&time,sun_point,lres)
    print(max(var1(nc,:,:)))
    print(min(var1(nc,:,:)))
end do
;    resp@txString      = "Year" + ny
    gsn_panel(wks,plot,(/ncase,1/),resp)
end

