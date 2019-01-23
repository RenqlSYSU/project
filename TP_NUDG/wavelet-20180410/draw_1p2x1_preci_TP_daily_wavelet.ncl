;*********************************************************************************
;filename:      draw_1p2x1_preci_TP_daily_wavelet.ncl
;input_data:    CTRL-preci_TP_JJAdaily_wavelet.nc
;                                            by Ql_Ren
;                                           2018/04/10
;*********************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.015
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

area =(/"STP","CTP","NTP","TP","SEA","NEA","EA"/)  ;S mean south,N mean North,EA mean East Asia
lats = (/25  ,30   ,35   ,25  , 22  , 30  , 20 /) ;-20
latn = (/30  ,35   ,40   ,40  , 30  , 45  , 45 /)  ;60
lonl = (/70  ,70   ,70   ,70  ,107  ,107  ,107 /)  ;20
lonr = (/105 ,105  ,105  ,105 ,130  ,130  ,122 /);130 ;220
lev  = (/850,500,200/)
draw_var   = (/"U","V","T","preci"/)
unit       = (/"m~S~2~N~/s~S~2~N~","m~S~2~N~/s~S~2~N~","k~S~2~N~","mm~S~2~N~/day~S~2~N~"/)
var_name   = (/"wave"/)
na = 3
nv = 2
nl = 2
isnormal = 0 ;1 mean normalized

if (isnormal.eq.1) then
cnlevels =  (/0.002,0.005 ,0.009,0.015,0.02 , \
              0.03 ,0.05  ,0.08 ,0.12 ,0.18 , \
              0.25 ,0.35  , 0.5 , 0.7 ,1    /)
else
cnlevels  = (/0.01,0.03,0.06,0.1 ,0.15, \
              0.2 ,0.3 ,0.4 ,0.6 ,0.8 , \
              1   ,1.2 ,1.5 ,2   ,2.5 , 3 , 4/)
end if

path    = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/wavelet-20180410/mdata/" + \
            (/"CTRL","NUDG6h","NUDG24h"/) + "_clim_daily_wavelet_" + area(na) + lev(nl) + draw_var(nv) + ".nc"
    case = (/"(a) CTRL", "(b) NG24h", "(c) NG6h"/)  
fig_out  = "pdf" ;"eps";"x11"

if(isnormal.eq.1) then
    fig_name   = path + "/wavelet-20180410/fig/clim_daily_wavelet_nobias-stand_" + area(na) + lev(nl) + draw_var(nv)
else
    fig_name   = path + "/wavelet-20180410/fig/clim_daily_wavelet_nobias_" + area(na) + lev(nl) + draw_var(nv)
end if

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
npart  = getfilevardimsizes(f,"part")
njtot  = getfilevardimsizes(f,"jtot")
ntime  = getfilevardimsizes(f,"time")
power  = new((/ncase,njtot,ntime/),float)
sig    = new((/ncase,njtot,ntime/),float)
plot   = new((/ncase/),graphic)
plot2  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
do nf = 0,ncase-1,1
    f   = addfile(filename(nf),"r")
    var = f->$var_name$       ;(ncase,neval,nlat,nlon)
    power(nf,:,:) = onedtond(var@power,(/njtot,ntime/))   
      sig(nf,:,:) = power(nf,:,:)/conform(power(nf,:,:),var@signif,0) ;compute significance ( >= 1 is significant)
    power(nf,:,:) = power(nf,:,:)/conform(power(nf,:,:),var@scale ,0)   ;no bias, referenced by Liu et.al(2007)
    if(isnormal.eq.1) then
        power(nf,:,:) = power(nf,:,:)/var@stdev/var@stdev
    end if
end do

power!1         = "period"
power&period    = var@period
power!2         = "time"
power&time      = var&time
copy_VarMeta(power,sig)
sig = mask(sig,sig.ge.1.0,True)
printVarSummary(power)
printVarSummary(sig)
printVarSummary(var)

exponent = ispan(2,8,2)
ylabel   = exp( exponent*log(2.0) )
print(ylabel)

;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@vpHeightF        = 0.2                   ; height and width of plot
    res@vpWidthF         = 0.9

;------------set the form of the coordinate-----------------------
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
    res@tmXTMajorOutwardLengthF = 0.01
    res@tmXBMajorOutwardLengthF = 0.01
    res@tmYLMajorOutwardLengthF = 0.01
    res@tmYRMajorOutwardLengthF = 0.01
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    res@tmXTMinorOutwardLengthF = 0.01
    res@tmXBMinorOutwardLengthF = 0.01
    res@tmYLMinorOutwardLengthF = 0.01
    res@tmYRMinorOutwardLengthF = 0.01

;---------set the attribute of the x y axis------------
    ;res@tiXAxisOn       = False
    res@tiXAxisString   = "Day"
    
    res@tiYAxisOn       = True
    res@tiYAxisFont     = font_index
    res@tiYAxisFontHeightF = subtitle_font 
    res@tiYAxisString   = "period(days)"
    
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = ispan(15,ntime,30)
    res@tmXBLabels      = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
    
    res@trYReverse          = True   ; reverse y-axis
    res@tmYLMode            = "Explicit"
    res@tmYLValues          = ylabel 
    res@tmYLLabels          = ylabel
    res@tmLabelAutoStride   = True

    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "black"
    res@tmYMajorGridThicknessF         = 4.0

    res_sig                  = True
    res_sig = res
;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
    res@cnFillMode          = "RasterFill"          ; turn on raster mode
    res@cnRasterSmoothingOn = True                  ; turn on raster smoothing

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    res@cnLevels             = cnlevels

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"   ;"Vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-2
    
;----------Set res sig-----------
    res_sig@tiXAxisOn       = False
    res_sig@tiYAxisOn       = False
    res_sig@tmLabelAutoStride   = True

    res_sig@cnFillOn         = True;False       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "grey0"
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    
    res_sig@gsnCenterString      = ""
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;------------set the attribute of the cone of influence---------------------
    rescoi              = True
    rescoi@gsFillIndex  = 1 

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString   = draw_var(nv) + " (" + unit(nv) +")"
    res@gsnCenterString  = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)" 

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

do nc = 0,ncase-1,1
    res@gsnLeftString    = case(nc)
    plot(nc)  = gsn_csm_contour(wks,power(nc,:,:),res)
    plot2(nc) = gsn_csm_contour(wks,sig(nc,:,:)  ,res_sig)
    overlay(plot(nc),plot2(nc))
    plot(nc) = ShadeCOI(wks, plot(nc), var, var&time, False)
    print(area(na))
    print(max(power(nc,:,:)))
    print(min(power(nc,:,:)))
end do

gsn_panel(wks,plot,(/ncase,1/),resp)
end

