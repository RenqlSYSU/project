;***************************************************************************
;filename:      draw_3p3x2_UVT_JJA_annual_ts.ncl
;plot the JJA times of every year's annual cycle(id=151~242) of TP average UVT(27 curves) 
;then overlaying the climate annual cycle(black bold curve)
;and one panel has six plots.
;                                            by Ql_Ren
;                                           2017/03/29
;******************************************************************************
begin
draw_var = (/"U"  ,"V"  ,"T"       ,"OMEGA"          ,"Z3"/)
unit     = (/"m/s","m/s","~S~o~N~C","10~S~-2~N~ Pa/s","gpm"/)   ;this can be used as right string 
var_name = (/"series","ave"/)
case     = (/" CTRL"," NG6h"/);this can be used as left string
pre_case = (/"(a","(b"/)
years    = ispan(1979,2005,1)   ;the calculate time

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/concept_map/mdata/" + (/"CTRL","NUDG6h"/) + "-" ;+ draw_var(nv) + "_TP_annual_series.nc"
fig_name   = path + "/concept_map/fig/5p3x2_UVTWZ_JJAts"
fig_out    ="pdf";"x11"

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.01
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 5.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(draw_var)
nlev  = dimsizes(lev)
nday  = 92
nday_s=151
nday_e=242
var   = new((/ncase,nvar,nyear,nday,nlev/),float)
ave   = new((/ncase,nvar,      nday,nlev/),float)
plot  = new((/ncase*nlev/),graphic)
plot2 = new((/ncase*nlev/),graphic)
plot3 = new((/ncase*nlev/),graphic)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nc = 0,ncase-1,1
do nv = 0,nvar-1,1
    f    = addfile(filename(nc)+draw_var(nv)+"_TP_annual_series.nc","r")
    var(nc,nv,:,:,:)  = f->$var_name(0)$(:,nday_s:nday_e,:)
    ave(nc,nv,:,:)    = f->$var_name(1)$(  nday_s:nday_e,:)
end do   ;the loop of the three variables 
end do   ;the loop of the three cases

var(:,2,:,:,:) = var(:,2,:,:,:) - 273.15   ;convert unit from k to c
ave(:,2,:,:)   = ave(:,2,:,:)   - 273.15
var(:,3,:,:,:) = var(:,3,:,:,:)*100   ;convert unit from k to c
ave(:,3,:,:)   = ave(:,3,:,:)*100

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
    ;resxy@tmYLAutoPrecision = False
    ;resxy@tmYLPrecision     = 2 
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    ;resxy@tmXBLabelFontHeightF     = tick_font-0.01
    ;resxy@tmXBLabelFontThicknessF  = tick_thick-1.0 
    ;resxy@tmYRLabelFontHeightF     = tick_font
    ;resxy@tmYRLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 90
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(0,95,5)
    resxy@tmXBLabels          = ispan(0,95,5)
    ;resxy@tmXBMode            = "Automatic"
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
    resxy@gsnRightStringFontHeightF  = subtitle_font
    resxy@gsnLeftStringFontHeightF   = subtitle_font
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

do nv =0, nvar-1,1
do nl =0, nlev-1,1
do nc = 0, ncase-1,1
    nn = nl +1
    np = ncase*nl + nc
    resxy@gsnLeftString   = pre_case(nc) + nn + case(nc)
    resxy@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"  
    resxy@trYMinF = floor(min(var(:,nv,:,:,nl)))
    resxy@trYMaxF =  ceil(max(var(:,nv,:,:,nl)))
    plot(np)  = gsn_csm_xy(wks,ispan(0,91,1),(/var(nc,nv,:,:,nl)/),resxy)
    plot2(np) = gsn_csm_xy(wks,ispan(0,91,1),(/ave(0 ,nv,  :,nl)/),resxy2)
    plot3(np) = gsn_csm_xy(wks,ispan(0,91,1),(/ave(nc,nv,  :,nl)/),resxy3)
    overlay(plot(np),plot2(np))                             
    overlay(plot(np),plot3(np))
end do
end do
    
   gsn_panel(wks,plot,(/nlev,ncase/),resp)
end do
end

