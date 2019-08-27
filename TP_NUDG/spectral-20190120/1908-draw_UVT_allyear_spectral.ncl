;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql"
var_name = (/"U","V","T"/)
case_option = 1 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
else
    ;remove TP experiment, sample number is ten years 
    case1   = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    case    = (/"CTRL","noTP"/)
    filein  = path + "/model/" + case1 + "/pro/" + case1 + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)

nv = 2
nl = 2
;nc = 0

;fig_name  = path + "/project/TP_NUDG/spectral-20190120/fig/" + case(nc) + "_" + var_name(nv) + "_" + lev(nl)
fig_name  = path + "/project/TP_NUDG/spectral-20190120/fig/1p3x3_" + var_name(nv)
;fig_name  = path + "/project/TP_NUDG/spectral-20190120/fig/3p3x3_UVT"
;fig_title = case(nc) + " " + var_name(nv) + " " + lev(nl)

f   = addfile(filein(0)+years(0)+".daily."+var_name(nv)+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime       = dimsizes(select_time)
ncase = dimsizes(case)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
alltime = ntime*nyear

;define array to calc the average of 27 years
var  = new((/nyear,ntime/),float)
plot = new((/ncase*nlev/),graphic)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

wks  = gsn_open_wks("pdf",fig_name)             ; send graphics to PNG file
np = 0
do nl = 0, nlev-1 ,1
do nc = 0, ncase-1 ,1
years   = ispan(21,30,1)   ;the calculate time
if((case_option.eq.1).and.(nc.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var1 = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})
    var(nt,:) = wgt_areaave_Wrap(var1,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
end do
;var = var - conform(var,dim_avg_n(var,0),1)  ;remove annual cycle

do nt = 0, nyear-1 ,1
; calculate spectrum
  sdof = specx_anal(ndtooned(var),d,sm,pct) ;returns the degrees of freedom as a scalar.

; return 4 curves to be plotted, first is the input spectrum, second is Markov "Red Noise" spectrum
; the last two numbers are the confidence levels
  splt = specx_ci (sdof, 0.05, 0.95) ;return splt(4,nyear*ntime/2)
end do
printVarSummary(sdof)
printVarSummary(splt)        ; splt(4,:)

frq = sdof@frq          ; convenience
p   = 1/frq             ; *highly non-linear*
p!0 = "frq"
p&frq = frq
p@long_name = "period"
p@units     = "day"
p_draw = (/60,20,10,7,5/)

i=0
j=1000
k=20
 print("====")
 print(frq(i:j:k)+"   "+p(i:j:k)+"  "+splt(0,i:j:k) )             ; print 
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
   mp_thick       = 2.0
   font_index     = 22

    resxy     = True
    ;resxy@vpWidthF  = 0.45
    ;resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = (/0,0,1,3/)
    resxy@xyLineThicknesses = (/2.0,3.0,3.0,3.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","red","red","red"/)
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "Period (day)"          ; xaxis
    resxy@tiYAxisString = "Variance/freq"           ; yaxis
    ;resxy@trYLog  = True
    resxy@trYMinF = 0
    resxy@trYMaxF = 500
    ;resxy@trXLog  = True
    ;resxy@trXMinF = 1
    ;resxy@trXMaxF = 10000
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 5
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 90
    ;resxy@tmXBMode            = "Explicit"
    ;resxy@tmXBValues          = 1.000/p_draw
    ;resxy@tmXBLabels          = p_draw
    
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
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = case(nc)
    resxy@gsnRightString   = var_name(nv) + " " + lev(nl) + "hPa" 

   ;ip   = ind(frq.le.0.2)    ; all indices for "short" periods
   ip   = ind(frq.ge.0.01)    ; all indices for "short" periods
   plot(np) = gsn_csm_xy(wks,p(ip), splt(:,ip),resxy)     
   ;plot(np) = gsn_csm_xy(wks,p, splt,resxy)     
   ;plot = gsn_csm_xy(wks, frq(ip), splt(:,ip), resxy)     
   ;plot = gsn_csm_xy(wks, frq, splt, resxy)     
   np = np + 1
   print(max(splt(:,ip)))
end do 
end do

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/nlev,ncase/),resp)

;draw(wks)
;frame(wks)
end
