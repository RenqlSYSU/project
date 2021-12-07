;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql/project/2021variability"
filein  = path + "/data/CN05.1-eof_ts_68_preci_5run.nc"
fig_name = path + "/fig/CN05.1-eof_ts_68_preci_5run"
pre_case = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out  ="pdf";"x11"
case2 = (/"OBS"   ,"OBS"  ,"CTL" ,"NT"   ,"NM"    ,"NMT" /) 
dcase = (/1,2/)
years = ispan(1979,2005,1)   ;the calculate time
nyear = dimsizes(years)
ntime = 92

;---------------------------------------------------------------
;read the data, then calculate the yearly power spectra 
;--------------------------------------------------------------
f    = addfile(filein,"r")
var  = f->eof_ts(:,0:1,:)

do nc = 0,1,1
    r = escorc(var(nc,0,:),var(nc,1,:))
    prob = rtest(r,toint(nyear*ntime/5),0)
    print(case2(dcase(nc))+"  r: "+r+"  prob: "+prob)
end do

var := dim_avg_n(reshape(var,(/2,2,nyear,ntime/)),2)

plot = new((/2/),graphic)

wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.013
   title_font     = 0.01
   label_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 2.0
   font_index     = 22

    resxy     = True
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    resxy@tiXAxisString = "day"          ; xaxis
    resxy@tiYAxisString = "PC"
    resxy@tiXAxisFontHeightF = tick_font
    resxy@tiYAxisFontHeightF = tick_font
    resxy@tiXAxisFont = font_index
    resxy@tiYAxisFont = font_index
    resxy@trYMinF = min(var)
    resxy@trYMaxF = max(var)
    resxy@trXMinF = 0 
    resxy@trXMaxF = ntime 
    resxy@xyComputeXMin = True
    ;resxy@tmXBMode    = "Explicit"
    ;resxy@tmXBValues  = ispan(0,nyear*ntime-1,ntime) 
    ;resxy@tmXBLabels  = years
    ;resxy@tmXBValues  = ispan(0,ntime-1,ntime) 
    ;resxy@tmXBLabels  = ispan(0,ntime-1,1)
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case2(dcase) 
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@lgLabelFont            = font_index
    resxy@pmLegendWidthF   = 0.1
    resxy@pmLegendHeightF  = 0.05
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = - 0.28   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.28    ; move ref vector along X-axis
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font 
    resxy@tmXBLabelFontThicknessF  = tick_thick
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

    ;resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    ;resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    ;resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    ;resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    ;resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    ;resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont        = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa"
    resxy@gsnRightString   = "" 

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    resxy@xyDashPatterns    = (/0,0/)
    resxy@xyLineThicknesses = (/2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"red","blue"/)
    
    np = 0
    do nc = 0, 1, 1
        resxy@gsnLeftString   = pre_case(np) + ") PC" + (np+1)
        plot(np) = gsn_csm_xy(wks,ispan(0,ntime-1,1), var(:,nc,:), resxy)     
        ;plot(np) = gsn_csm_xy(wks,ispan(0,nyear*ntime-1,1), var(:,nc,:), resxy)     
        np = np + 1
    end do
    gsn_panel(wks,plot,(/2,1/),resp)
end
