;***************************************************************************
;filename:      draw_1p3x3_UVT_annual_series_MethodologyPaper.ncl
;plot the every year's annual cycle(356) of TP average UVT(27 curves) 
;then overlaying the climate annual cycle(black bold curve)
;and one panel has nine plots
;this fig is used to the paper.
;                                            by Ql_Ren
;                                           2017/04/11
;******************************************************************************
begin
draw_var = (/"corr"/)
var_name = (/"corr"/)
case1    = (/"CTRL","NUDG6h"/)
case2    = (/"CTRL","NG6h-CTRL"/)
month_s = 3
month_e = 7
mxlag    = 10
nc = 1

;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
filein = path + "/correlation-20191018/mdata/" + case1(nc) + "-daily" + mxlag + "delay_corri_" + month_s + "-" + month_e + ".nc"
fig_name = path + "/correlation-20191018/fig/" + case1(nc) + "-daily" + mxlag + "delay_corri_" + month_s + "-" + month_e
fig_title= "Daily Delay Correlation" 
fig_out  = "pdf" ;"eps" ;"x,11"
pre_case = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)

varname = (/"preci","QTE850","EKE850","EKE200","q1_850","q1_200","div850","div200"/)
nvar2 = dimsizes(varname)
nvar3 = (nvar2-1)*nvar2/2
label = new((/nvar3/),string)
npanel = 7
nv3 = 0
do nv1 = 0,nvar2-2, 1
do nv2 = nv1+1, nvar2-1, 1
    label(nv3) = varname(nv2) + "&" + varname(nv1)
    nv3 = nv3 + 1
end do
end do

area = (/"SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 22  , 30  /)
latn = (/ 30  , 45  /)
lonl = (/107  ,115  /)
lonr = (/130  ,130  /)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.013
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

;==========================================================
;read data
;============================================================
f = addfile(filein,"r")
var   = f->$var_name$
print(var(0,2,:))
narea = dimsizes(area)
plot  = new((/narea/),graphic)

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    resxy     = True
    resxy@vpWidthF  = 0.4
    resxy@vpHeightF = 0.3
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = (/0,0,0,0/) ;,0,0
    resxy@xyLineThicknesses = (/4.0, 4.0, 4.0, 4.0/) ;, 4.0, 4.0  ; make second line thicker
    resxy@xyLineColors      = (/"black","red","blue","grey"/); ,"darkorange3","chartreuse4"
    resxy@tiXAxisOn = True
    resxy@tiXAxisString = "lag(days)"
    resxy@tiYAxisOn = False
    resxy@trYMinF = -0.5 ;floor(min(var)); -0.4;
    resxy@trYMaxF =  0.5 ;ceil(max(var)) ; 0.2 ;
    
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
    ;resxy@tmXBTickSpacingF   = 30
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 360
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(-10,10,5)
    resxy@tmXBLabels          = ispan(-10,10,5)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    resxy@gsnYRefLine           = 0.0
    resxy@gsnYRefLineColor      = "black"
    resxy@gsnYRefLineThicknessF = 5.0
    resxy@gsnYRefLineDashPattern= 0

    ;legend
    resxy@pmLegendDisplayMode = "NoCreate";"Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    ;resxy@xyExplicitLegendLabels = label
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.09
    resxy@pmLegendHeightF  = 0.3
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.05   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.3    ; move ref vector along X-axis
    
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF   = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = "" 
    resxy@gsnRightString   = ""

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

do np = 0,npanel-1,1
    nvs = np*(nvar3/npanel)
    nve = (np+1)*(nvar3/npanel)-1
    resxy@xyExplicitLegendLabels = label(nvs:nve)
    do na = 0, narea-1,1
        resxy@gsnLeftString   = pre_case(na) + ") " + case2(nc)
        resxy@gsnRightString  = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"  
        if(na.eq.1) then
        resxy@pmLegendDisplayMode = "Always"
        end if
        plot(na)  = gsn_csm_xy(wks,var&delay,(/var(na,nvs:nve,:)/),resxy)
    end do
    
   gsn_panel(wks,plot,(/1,2/),resp)
end do
end


