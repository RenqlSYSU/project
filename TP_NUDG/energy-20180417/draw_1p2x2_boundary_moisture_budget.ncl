;*********************************************************************
;filename:      draw_1p2x1_boundary_moisture_budget.ncl
;input_data:    CTRL-Clim_JJAave_qtran_int.nc
;draw four plots(upside,bottom,left,right), every plot have two curves(CTRL,NG6h),
;only draw climate transport
;                                                     by Ql Ren
;                                                    2018/04/24
;*********************************************************************************
begin
draw_var   = (/"climate"/) 
var_name_u = (/"qu"/)
var_name_v = (/"qv"/)
unit     = "kg/(s*m)" ;"K"   ;this can be used as right string 

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h"/) + "-Clim_JJAave_qtran_int.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/1p2x2_boundary_moisture_budget"
    case   = (/"CTRL","NG6h"/)
plot_right = (/"upside","bottom","left","right"/)
pre_case   = (/"(a)","(b)","(c)","(d)"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = 25  ;-15 ;-20
latn = 40  ;55  ;60
lonl = 70  ;30  ;20
lonr = 105 ;150 ;220
siglvl = 0.05

f      = addfile(filename(0),"r")
vars  := f->$var_name_u(0)$({lats:latn},{lonl:lonr})
nplot  = dimsizes(plot_right)
ncase  = dimsizes(case)
nlat   = dimsizes(vars&lat)
nlon   = dimsizes(vars&lon)
bottom = new((/ncase,nlon/),float)
upside = new((/ncase,nlon/),float)
left   = new((/ncase,nlat/),float)
right  = new((/ncase,nlat/),float)

;the number of the panel plot
    plot     = new((/nplot/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.025
   subtitle_font  = 0.03
   title_font     = 0.03
   label_font     = 0.025
   vcref_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   vc_thick       = 1.5
   font_index     = 22

;================================================================
;read data
;====================================================================
do nf = 0,ncase-1, 1
    f = addfile(filename(nf),"r")
    bottom(nf,:) = f->$var_name_v$({lats},{lonl:lonr})
    upside(nf,:) = f->$var_name_v$({latn},{lonl:lonr})
      left(nf,:) = f->$var_name_u$({lats:latn},{lonl})
     right(nf,:) = f->$var_name_u$({lats:latn},{lonr})
end do

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.4
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
   
    resxy@xyDashPatterns    = (/0 ,0 /)
    resxy@xyLineThicknesses = (/5.0, 5.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"red","blue"/)
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
    resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 360
    ;resxy@tmXBMode            = "Explicit"
    ;resxy@tmXBValues          = ispan(0,360,30)
    ;resxy@tmXBLabels          = ispan(0,360,30)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnCenterString      = ""
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF  = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 

    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case 
    resxy@lgLabelFontColor       = "black"
    ;resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@lgLabelFont            = font_index
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.08
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = -0.4   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.4   ; move ref vector along X-axis

do np = 0, nplot-1,1
    resxy@gsnLeftString    = pre_case(np) +" " + plot_right(np)
    resxy@gsnRightString   = unit
    if(np.le.1) then
    resxy@tmXBTickSpacingF   = 5
    if(np.eq.0) then
        plot(np)  = gsn_csm_xy(wks,upside&lon,(/upside/),resxy)
    else
        plot(np)  = gsn_csm_xy(wks,bottom&lon,(/bottom/),resxy)
    end if
    else
    resxy@tmXBTickSpacingF   = 3
    resxy@trXMaxF   = 40
    resxy@trXMinF   = 25
    if(np.eq.2) then
        plot(np)  = gsn_csm_xy(wks,left&lat,(/left/),resxy)
    else
        plot(np)  = gsn_csm_xy(wks,right&lat,(/right/),resxy)
    end if
    end if
end do
    
    resp  = True    
    resp@gsnMaximize   = True
   gsn_panel(wks,plot,(/2,2/),resp)
end
