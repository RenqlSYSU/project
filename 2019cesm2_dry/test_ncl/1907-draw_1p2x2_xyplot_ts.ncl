;*********************************************************************
;draw temperature (color) and pressure (contour) at every time and level
;                                                     by Ql Ren
;                                                    2019/06/18
;*********************************************************************************
begin
draw_var = (/"T","QRS"/)  ;this is the variable used to draw, can used as center string
unit     = (/"k","K/S"/)    ;this can be used as right string 
lev = (/5,15,29/) 
nl  = 2
nv  = 1

;area: east Asia
lats = -90 ;20 ;-20
latn =  90 ;55 ;60
lonl = 0   ;30 ;20
lonr = 360 ;150;220
lat_spacing = 30
lon_spacing = 60
lon  = (/0,90,180,270/)
nlon = dimsizes(lon)

;the path of the file
path   = "/home/ys17-19/renql"
case1  = (/"eul_HS94_2"/)
case2  = (/"eul_HS94_2"/)
nf = 0

filename = path + "/model/" + case1(nf) + "/" + case1(nf) + ".cam.h0.*"
fig_name = path + "/project/2019cesm2_dry/test_ncl/fig/1p2x2_" + case1(nf) + "_xyplot_"  + lev(nl) + draw_var(nv)
;fig_title  = case2(nf) + " lev = " + lev(nl)
;fig_title  = case2(nf) + " lev = " + lev(nl)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

;the number of the panel plot
    nplot    = 4
    plot     = new((/nplot/),graphic)
    plot2    = new((/nplot/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.5
   font_index     = 22

;================================================================
;read data
;====================================================================
files       = systemfunc("ls " + filename) 
f           = addfiles(files,"r") ;(9:31)
ListSetType(f,"cat")
times       = f[:]->time
ntimes = dimsizes(times)         ; number of times in the file
printVarSummary(times)

var = f[:]->$draw_var(nv)$(:,lev(nl),:,{lon})
ave = dim_avg_n_Wrap(var,0)
printVarSummary(var)
printVarSummary(ave)

lat = var&lat*3.1415926/180.0
psurf_ref = f[:]->P0(0) 
p = f[:]->lev(lev(nl))*100 ;unit is Pa
print(psurf_ref)
print(p)
cappa = 2.0 / 7.0
trefc = 315 - (60 * sin(lat))
trefa = (trefc - (10*(cos(lat)^2)*log(p/psurf_ref)))*(p/psurf_ref)^cappa
trefa = where(trefa.lt.200,200,trefa)

wks = gsn_open_wks(fig_out,fig_name)

;==============================================================
;draw figure at every time  
;==================================================================
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    resxy     = True
    resxy@vpWidthF  = 0.45
    resxy@vpHeightF = 0.25
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = ispan(2,64,2)
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
    ;resxy@tmXBMode            = "Explicit"
    ;resxy@tmXBValues          = ispan(0,360,30)
    ;resxy@tmXBLabels          = ispan(0,360,30)
    
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

    ;resxy@tiXAxisOn = True
    ;resxy@tiXAxisFont        = font_index
    ;resxy@tiXAxisFontHeightF = tick_font
    ;resxy@tiXAxisString      = "day"

do nlo = 0, nlon-1,1
    resxy@gsnLeftString    = lon(nlo)
    resxy@gsnRightString   = draw_var(nv) + " (" + unit(nv) + ")"
    resxy@trYMinF = floor(min(var(:,:,:)))
    resxy@trYMaxF =  ceil(max(var(:,:,:)))
    plot(nlo)  = gsn_csm_xy(wks,var&lat,(/var(:,:,nlo)/),resxy)
    if(nv.eq.0) then 
    plot2(nlo) = gsn_csm_xy(wks,var&lat,trefa,resxy2)
    else
    plot2(nlo) = gsn_csm_xy(wks,var&lat,(/ave(:,nlo)/),resxy2)
    end if
    overlay(plot(nlo),plot2(nlo))
end do
fig_title  = case2(nf) + " lev = " + p + "Pa"
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@txString      = fig_title
   gsn_panel(wks,plot,(/2,2/),resp)
end

