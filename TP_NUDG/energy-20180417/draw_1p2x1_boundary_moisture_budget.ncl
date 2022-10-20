;*********************************************************************
;filename:      draw_1p2x1_boundary_moisture_budget.ncl
;input_data:    CTRL-Clim_JJAave_qtran_int.nc 
;                                                     by Ql Ren
;                                                    2018/04/24
;*********************************************************************************
begin
draw_var   = (/"climate","steady","transient"/) 
var_name_u = (/"qu","ste_qu","edd_qu"/)
var_name_v = (/"qv","ste_qv","edd_qv"/)
unit     = "kg/(s*m)" ;"K"   ;this can be used as right string 

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h"/) + "-Clim_JJAave_qtran_int.nc"
fig_name = path + "/energy-20180417/fig/1p2x1_boundary_moisture_budget"
txt_name = path + "/energy-20180417/TP_boundary_moisture_budget.txt"
    case   = (/"CTRL","NG6h"/)
plot_right = (/"bottom_upside","left_right"/)
pre_case   = (/"(a","(b"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = 25  ;-15 ;-20
latn = 40  ;55  ;60
lonl = 70  ;30  ;20
lonr = 105 ;150 ;220
siglvl = 0.05

f      = addfile(filename(0),"r")
vars  := f->$var_name_u(0)$({lats:latn},{lonl:lonr})
nvar   = dimsizes(draw_var)
ncase  = dimsizes(case)
nlat   = dimsizes(vars&lat)
nlon   = dimsizes(vars&lon)
bottom = new((/ncase,nvar,nlon/),float)
upside = new((/ncase,nvar,nlon/),float)
left   = new((/ncase,nvar,nlat/),float)
right  = new((/ncase,nvar,nlat/),float)
up_bot  = new((/ncase*nvar*2,nlon/),float)
lef_rit = new((/ncase*nvar*2,nlat/),float)
label_string = new((/2,ncase*nvar*2/),string)

;the number of the panel plot
    plot     = new((/2/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.02
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
do nv = 0,nvar-1,1
    f = addfile(filename(nf),"r")
    bottom(nf,nv,:) = f->$var_name_v(nv)$({lats},{lonl:lonr})
    upside(nf,nv,:) = f->$var_name_v(nv)$({latn},{lonl:lonr})
      left(nf,nv,:) = f->$var_name_u(nv)$({lats:latn},{lonl})
     right(nf,nv,:) = f->$var_name_u(nv)$({lats:latn},{lonr})
end do
end do

nn = 0
do nf = 0,ncase-1, 1
do nv = 0,nvar-1,1
    label_string(0,nn) = case(nf) + "_bottom_" + draw_var(nv)
    label_string(1,nn) = case(nf) + "_left_" + draw_var(nv)
    up_bot(nn,:)  = bottom(nf,nv,:)
    lef_rit(nn,:) = left(nf,nv,:)
    nn = nn +1
    label_string(0,nn) = case(nf) + "_upside_" + draw_var(nv)
    label_string(1,nn) = case(nf) + "_right_" + draw_var(nv)
    up_bot(nn,:)  = upside(nf,nv,:)
    lef_rit(nn,:) = right(nf,nv,:)
    nn = nn +1
end do
end do

dlon = 1.25 ;unit is angle
dlat = 0.9  ;unit is angle
earth_radius = 6378388.0 ;unit is m
angle2radian = 3.14159265/180  ;arc length divided by radius is radian,such as pi
dxs = earth_radius*cos(angle2radian*lats)*(dlon*angle2radian)  ;unit is m
dxn = earth_radius*cos(angle2radian*latn)*(dlon*angle2radian)
dy  = earth_radius*(dlat*angle2radian)

bave = dim_sum_n(bottom*dxs,2)/1000000   ;(ncase,nvar), unit is kg/s
uave = dim_sum_n(-upside*dxn,2)/1000000  ;positive means the moisture move in
lave = dim_sum_n(left*dy,2)/1000000      
rave = dim_sum_n(-right*dy,2)/1000000    
printVarSummary(bave)
printVarSummary(bottom)
print(dy)
print(dxs)

;-----write boundary moisture anlysis result to txt------------------
stri = (/"Climate","steady","transient"/)
write_table(txt_name,"w",[/"bottom","upside","left","right","unit is 10^6 kg/s"/],"%s   %s  %s   %s   %s")
do nc = 0,ncase-1,1
    ave = [/stri,bave(nc,:),uave(nc,:),lave(nc,:),rave(nc,:)/]
    write_table(txt_name,"a",ave,"%s %16.4f %16.4f %16.4f %16.4f")
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
   
    resxy@xyMarkLineMode    = "MarkLines"
    resxy@xyMarkers         = (/16,4 ,16,4 ,16,4 , 16,4 ,16,4 ,16,4/)
    resxy@xyDashPatterns    = (/0 ,0 ,0 ,0 ,0 ,0 , 5 ,5 ,5 ,5 ,5 ,5/)
    resxy@xyLineThicknesses = 5.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","black","red","red","blue","blue","black","black","red","red","blue","blue"/)
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    
    resxy@tmBorderThicknessF  = mp_thick
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
    ;resxy@xyExplicitLegendLabels = label_string 
    resxy@lgLabelFontColor       = "black"
    ;resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@lgLabelFont            = font_index
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.4
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.7   ; move ref vector along X-axis

do nc = 0, ncase-1,1
    resxy@gsnLeftString    = pre_case(nc) + ") " + plot_right(nc)
    resxy@gsnRightString   = unit
    resxy@xyExplicitLegendLabels = label_string(nc,:)
    if(nc.eq.0) then
    resxy@tmXBTickSpacingF   = 5
    plot(nc)  = gsn_csm_xy(wks,bottom&lon,(/up_bot/),resxy)
    else
    resxy@tmXBTickSpacingF   = 3
    plot(nc)  = gsn_csm_xy(wks,left&lat,(/lef_rit/),resxy)
    end if
end do
    
    resp  = True    
    resp@gsnMaximize   = True
   gsn_panel(wks,plot,(/2,1/),resp)
end
