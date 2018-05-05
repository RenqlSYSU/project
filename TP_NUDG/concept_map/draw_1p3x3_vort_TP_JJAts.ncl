;***************************************************************************
;filename:      calc_clim_month_ave_vort.ncl
;output_file:   Clim_month_ave_vort.nc 
;read the 27 years variable(TLLL) of one case at East AsiaI(15S£­55N£¬30-150E) 
;and then calculate the month average to get an array of (nmonth,nlat,nlon),
;then calculate the season average (nseason,nlat,nlon).
;                                            by Ql_Ren
;                                           2018/01/18
;******************************************************************************
begin
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

path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fig_name = path + "/concept_map/fig/1p5x3_vort_SCNTP_SNEA_JJAts"
fig_out  = "eps";"pdf";"x11"
;fig_out  = "pdf";"x11"

case    = (/"CTRL","NG24h","NG6h"/)
var_name  = (/"U","V"/)
calc_name = (/"vort"/)
unit     = (/"10~S~-6~N~ s~S~-1~N~"/)   ;this can be used as right string 
pre_case = (/"(a","(b","(c","(d","(e"/)
years    = ispan(1979,2005,1)   ;the calculate time

area = (/"STP","CTP","NTP","SEA","NEA"/);"STP","CTP","NTP","SEA","NEA"/) "TP", ;S mean south,N mean North,EA mean East Asia
lats = (/25   ,30   ,35   , 25  , 35  /);25   ,30   ,35   , 22  , 30  /) 25  , ;-20
latn = (/30   ,35   ,40   , 35  , 45  /);30   ,35   ,40   , 30  , 45  /) 40  , ;60
lonl = (/70   ,80   ,70   ,107  ,107  /);70   ,70   ,70   ,107  ,107  /) 70  , ;20
lonr = (/105  ,100  ,105  ,120  ,130  /);105  ,105  ,105  ,130  ,130  /) 105 , ;130 ;220
lev = (/500/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
narea = dimsizes(area)
nday  = 90
nday_s=152
nday_e=241
var   = new((/ncase,narea,nyear,nday/),float)
ave   = new((/ncase,narea,      nday/),float)
plot  = new((/ncase*narea/),graphic)
plot2 = new((/ncase*narea/),graphic)
plot3 = new((/ncase*narea/),graphic)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nc = 0,ncase-1,1
do na = 0,narea-1,1
do nt = 0,nyear-1,1
    nv = 0
    f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_u := f->$var_name(nv)$(nday_s:nday_e,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})

    nv = 1
    f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_v := f->$var_name(nv)$(nday_s:nday_e,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})

    vort := uv2vr_cfd(var_u, var_v, var_v&lat, var_v&lon, 0) ;calculate the vorticity,get(nday,nlat,nlon)
    printVarSummary(vort)
    var(nc,na,nt,:) = wgt_areaave(vort,1.0,1.0,0)
end do
end do
end do
var = var*1000000
ave = dim_avg_n_Wrap(var,2)
printVarSummary(var)
printVarSummary(ave)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    resxy     = True
    resxy@vpWidthF  = 0.45
    resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = ispan(4,56,2)
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    
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
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 90
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(0,90,10)
    resxy@tmXBLabels          = ispan(0,90,10)
    resxy@trXMaxF      = 90
    resxy@trXMinF      = 0
    ;resxy@tmXBMode            = "Automatic"
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""

    resxy2 = resxy
    resxy2@xyMonoLineColor  = True
    resxy2@xyMonoLineThickness = True
    resxy2@xyLineThicknessF = 4.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy2@xyLineColor      = "black"

    resxy3 = resxy2
    resxy3@xyLineThicknessF = 6.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

do na =0, narea-1,1
do nc = 0, ncase-1,1
    nn = nc +1
    np = ncase*na + nc
    resxy@gsnLeftString   = pre_case(na) + nn + ")" + case(nc)
    resxy@gsnRightString  = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"  
    resxy@trYMinF = floor(min(var(:,na,:,:)))
    resxy@trYMaxF =  ceil(max(var(:,na,:,:)))
    plot(np)  = gsn_csm_xy(wks,ispan(0,89,1),(/var(nc,na,:,:)/),resxy)
    ;plot2(np) = gsn_csm_xy(wks,var&day,(/ave(0 ,na,  :)/),resxy2)
    plot3(np) = gsn_csm_xy(wks,ispan(0,89,1),(/ave(nc,na,  :)/),resxy3)
    ;overlay(plot(np),plot2(np))
    overlay(plot(np),plot3(np))
end do
end do
    
   gsn_panel(wks,plot,(/narea,ncase/),resp)
end

