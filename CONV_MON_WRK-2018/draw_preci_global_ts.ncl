;***************************************************************************
;filename:      draw_fsnt_flnt_preci_global_ts.ncl
;1 read data in (ncase,12nmonth*100years,nlat,nlon),FSNT-FLNT
;2 calc global weighted annual mean，then draw the ts plot
;                                            by Ql_Ren
;                                           2018/07/19
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
draw_var = (/"preci"/)  ;the radiation at top of atmosphere,this is the variable used to draw, can used as center string
unit     = (/"mm/day"/)    ;this can be used as right string 
var_name = (/"PRECC","PRECL"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)

;the path of the file
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fig_name = path + "/fig/" + draw_var + "_global_ts"
fig_out  = "pdf";"x11"

f    = addfile(filename(0) + var_name(0) + ".1-100.nc","r")
vars = f->$var_name(0)$(0,:,:)
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
ncase  = dimsizes(case)
nyear  = 100 
nmonth = 12
nvar   = dimsizes(var_name)
var   = new((/ncase,100*12,nlat,nlon/),float)
var1  = new((/ncase,nyear /),float)
;plot  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
    nv = 0
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    fsnt =  f->$var_name(nv)$ ;+ means downward to earth
    
    nv = 1
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    flnt =  f->$var_name(nv)$ ;+ means upward to space

    var(nc,:,:,:) = fsnt+flnt ;+ means downward to earth
end do
var = var*3600*24*1000 ;convert unit from m/s to mm/day
copy_VarMeta(vars,var(0,0,:,:))
printVarSummary(var)

pai  = 3.1415926
re   = 6371220.0
dlon = 2.5
lat  = vars&lat
dx = re*cos(lat*pai/180.0)*(dlon*pai/180.0)
wgty = dx
wgtx = 1.0
var_area = wgt_areaave(var,wgty,wgtx,0)

do ns = 0,nyear-1,1
    nms = ns*nmonth
    nme = ns*nmonth + 11 
    var1(:,ns) = dim_avg_n(var_area(:,nms:nme),1) ;calc the annual mean
end do
printVarSummary(var1)

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
   mp_thick       = 4.0
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.3
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      =(/"black","red","green","blue"/)
    resxy@tiXAxisOn = True
    resxy@tiXAxisString = "Year"
    resxy@tiYAxisOn = False
    resxy@trYMinF = 2.8;floor(min(var1))
    resxy@trYMaxF = 3.6; ceil(max(var1))
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font 
    resxy@tmXBLabelFontThicknessF  = tick_thick
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    resxy@tmXBMode           = "Manual"
    resxy@tmXBTickSpacingF   = 10
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 100
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    resxy@gsnYRefLine            = 0
    resxy@gsnYRefLineDashPatterns = 12 ; (/12,12,12/)   
    resxy@gsnYRefLineThicknesses  = 4 ; (/3.0,3.0,3.0/)
    resxy@gsnYRefLineColors       = "black" ; (/"blue","red","green"/)
                             
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.2
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = -0.3   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.35    ; move ref vector along X-axis
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font

    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = "global " + draw_var + " (" + unit + ")"
    resxy@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape"

;    resxy@gsnLeftString   = pre_case(nc) + ") " + case(nc)
    plot  = gsn_csm_xy(wks,ispan(1,100,1),(/var1/),resxy)
    
   gsn_panel(wks,plot,(/1,1/),resp)
end


