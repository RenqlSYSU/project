;***************************************************************************
;filename:      draw_1p1x2_intro_sketch_ts.ncl
;use the daily annual cycle of TP area average T in 850hPa 0f CTRL and NG6h 
;these two case to explain the mean of nudging in the intro of the paper.
;                                            by Ql_Ren
;                                           2018/04/11
;******************************************************************************
begin
case     = (/"(a) CTRL","(b) NUDG"/);this can be used as left string
years    = ispan(1979,2005,1)   ;the calculate time
draw_var = (/"U","V","T","OMEGA","Z3"/)
unit     = (/"m/s","m/s","~S~o~N~C","10~S~-2~N~ Pa/s","gpm"/)   ;this can be used as right string 
var_name = (/"series","ave"/)
nv = 2;3;1

;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/concept_map/mdata/" + (/"CTRL","NUDG6h"/) + "-" + draw_var(nv) + "_TP_annual_series.nc"
fig_name   = path + "/concept_map/fig/1p1x2_intro_sketch_ts" + draw_var(nv)
fig_out    ="pdf";"x11" ;"eps";

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)
nl   = 2

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
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
var         = new((/ncase,nyear,350/),float)
ave         = new((/ncase,      350/),float)
plot        = new((/ncase/),graphic)
plot2       = new((/ncase/),graphic)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nc = 0,ncase-1,1
    f    = addfile(filename(nc),"r")
    var(nc,:,:)  = f->$var_name(0)$(:,0:349,nl)
    ave(nc,:)    = f->$var_name(1)$(0:349,nl)
end do   ;the loop of the three cases

if (nv.eq.2) then
    var = var - 273.15   ;convert unit from k to c
    ave = ave - 273.15
end if

if (nv.eq.3) then
    var = var*100   ;convert unit from k to c
    ave = ave*100
end if

printVarSummary(var)
printVarSummary(ave)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    resxy     = True
    resxy@vpWidthF  = 0.6
    resxy@vpHeightF = 0.2
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = ispan(4,56,2)
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont   = font_index
    resxy@tmYLLabelFont   = font_index
    ;resxyxy@tmXBLabelFontHeightF     = tick_font
    ;resxyxy@tmYRLabelFontHeightF     = tick_font

    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    ;  resxy@tmXBMajorOutwardLengthF = -0.
    ;  resxy@tmYLMajorOutwardLengthF = -0.
    ;  resxy@tmYRMajorOutwardLengthF = -0.
    
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;  resxy@tmXBMinorOutwardLengthF = -0.
    ;  resxy@tmYLMinorOutwardLengthF = -0.
    ;  resxy@tmYRMinorOutwardLengthF = -0.

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 30
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 365
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(0,360,30)
    resxy@tmXBLabels          = ispan(0,360,30)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF  = subtitle_font 
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""

    resxy2 = resxy
    resxy2@xyMonoLineColor  = True
    resxy2@xyMonoLineThickness = True
    resxy2@xyLineThicknessF = 6.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy2@xyLineColor      = "black"

    resxy@tiXAxisOn = True
    resxy@tiXAxisFont        = font_index
    resxy@tiXAxisFontHeightF = tick_font
    resxy@tiXAxisString      = "day"

do nc = 0, ncase-1,1
    resxy@gsnLeftString    = case(nc)
    resxy@trYMinF = floor(min(var(:,:,:)))
    resxy@trYMaxF =  ceil(max(var(:,:,:)))
    plot(nc)  = gsn_csm_xy(wks,var&day,(/var(nc,:,:)/),resxy)
    plot2(nc) = gsn_csm_xy(wks,var&day,(/ave(0 ,  :)/),resxy2)
    overlay(plot(nc),plot2(nc))
end do
    
    resp  = True    
    resp@gsnMaximize   = True
   gsn_panel(wks,plot,(/1,ncase/),resp)
end

