;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h) and then 
;remove the annual cycle, then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql"
filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
var_name = (/"U","V","T"/)
case     = (/"CTRL","NG24h","NG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)

nv = 2
nl = 2
nc = 0

fig_name  = path + "/project/TP_NUDG/20190120-spectral/fig/" + case(nc) + "_" + var_name(nv) + "_" + lev(nl)
fig_title = case(nc) + " " + var_name(nv) + " " + lev(nl)

f   = addfile(filein(0)+years(0)+".daily."+var_name(nv)+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
ncase = dimsizes(case)
nyear = dimsizes(years)

;define array to calc the average of 27 years
var  = new((/nyear,ntime/),float)
splt = new((/4,toint(ntime/2)/),float)
splt = 0

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 5  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0.1;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

;do nc = 0, ncase-1 ,1
do nt = 0, nyear-1 ,1
f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
var1 = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})
var(nt,:) = wgt_areaave_Wrap(var1,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
end do
var = var - conform(var,dim_avg_n(var,0),1)  ;remove annual cycle

do nt = 0, nyear-1 ,1
; calculate spectrum
  sdof = specx_anal(var(nt,:),d,sm,pct) ;returns the degrees of freedom as a scalar.
; return 4 curves to be plotted, the last two numbers are the confidence levels
  splt_term = specx_ci (sdof, 0.05, 0.95)
  splt = splt + splt_term
end do
splt = splt/nyear
printVarSummary(sdof)
printVarSummary(splt_term)        ; splt(4,:)

frq = sdof@frq          ; convenience
p   = 1/frq             ; *highly non-linear*
p!0 = "frq"
p&frq = frq
p@long_name = "period"
p@units     = "day"
p_draw = (/60,20,10,7,5/)

 print("====")
 print(frq+"   "+p+"  "+splt(0,:) )             ; print 
 print("====")

;************************************************
; plotting
;************************************************
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

wks  = gsn_open_wks("png",fig_name)             ; send graphics to PNG file

    resxy     = True
    ;resxy@vpWidthF  = 0.45
    ;resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = (/0,0,1,3/)
    resxy@xyLineThicknesses = 3.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","red","red","red"/)
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "Period (day)"          ; xaxis
    resxy@tiYAxisString = "Variance/freq"           ; yaxis
    resxy@trYMinF = 0 
    resxy@trYMaxF = 20 
    
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
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 90
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = 1.000/p_draw
    resxy@tmXBLabels          = p_draw
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = "sm = " + sm
    resxy@gsnRightString   = fig_title

   ip   = ind(frq.le.0.2)                        ; all indices for "short" periods
   ;plot = gsn_csm_xy(wks,p(ip), splt(:,ip),res)     
   plot = gsn_csm_xy(wks, frq(ip), splt(:,ip), resxy)     

draw(wks)
frame(wks)
end
