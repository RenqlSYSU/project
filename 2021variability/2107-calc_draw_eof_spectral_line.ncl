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
var_name = (/"preci"/)
unit     = (/"mm/day"/)
nv = 0
ifrunave= 1 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
nave = 5 ;five days run ave
ca = 80.0 ;day
cb = 8.0
rmv = False;True ;remove annaul cycle when bandpass filter
allcase = False ;True ;if True, four case be draw in one figure
month_s = 6
month_e = 8

siglvl = 0.05
p_draw = (/120,90,60,40,30,20,15,10,7,5/)

case1   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
case2   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
years   = ispan(1979,2005,1)   ;the calculate time
if(ifrunave.eq.2)
if(rmv) then
tsfile  = path + "/data/"+case1(0)+"-eof_ts_"+month_s+month_e+"_"+var_name(nv)+"_"+cb+"-"+ca+"-rmv.nc"
else
tsfile  = path + "/data/"+case1(0)+"-eof_ts_"+month_s+month_e+"_"+var_name(nv)+"_"+cb+"-"+ca+".nc"
end if
fig_name = path + "/fig/" +case1(0)+"-eof_ts_"+month_s+month_e+"_"+cb+"-"+ca+"-spectral"
end if
if(ifrunave.eq.1) then 
tsfile  = path + "/data/"+case1(0)+"-eof_ts_"+month_s+month_e+"_"+var_name(nv)+"_"+nave+"run.nc"
fig_name = path + "/fig/" +case1(0)+"-eof_ts_"+month_s+month_e+"_"+nave+"run-spectral"
end if
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0.1;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

tsf = addfile(tsfile,"r")
eval = tsf->eval
pcvar= tsf->pcvar
ts = tsf->eof_ts ;ncase,neval,nyear*ntime
nyear = dimsizes(years)
dim   = dimsizes(ts)
ncase = dim(0)
neval = dim(1)
ntime = dim(2)

indx = 4 ;the result of specx_ci
nsample = toint(ntime/nyear)
nfrq = toint(nsample/2.0) 
splt = new((/ncase,indx,nfrq/),float)
eofts= reshape(ts,(/ncase,neval,nyear,nsample/))

wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file
do ne = 0, neval-2, 1
np = 0
nee = ne+1
do nc = 0, ncase-1, 1
    spcavg = new ((/nfrq/),typeof(eofts))
    spcavg = 0.0
    r1zsum = 0.0
    do nt = 0, nyear-1, 1
        ; calculate spectrum
        sdof := specx_anal(eofts(nc,ne,nt,:),d,sm,pct) ;returns the degrees of freedom as a scalar.
        spcavg = spcavg + sdof@spcx                ; sum spc of each segment
        r1     = sdof@xlag1                        ; extract segment lag-1
        r1zsum = r1zsum  + 0.5*log((1+r1)/(1-r1)) ; sum the Fischer Z
    end do 
    r1z  = r1zsum/nyear                 ; average r1z
    r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1); transform back, this is the mean r1
    spcavg  = spcavg/nyear              ; average spectrum
    df      = 2.0*nyear   ; deg of freedom
    df@spcx = spcavg      ; assign the mean spc
    df@frq  = sdof@frq
    df@xlag1= r1          ; assign mean lag-1
    splt(nc,:,:) = specx_ci(df, siglvl, 1.0-siglvl) 
end do
varia = splt*(1.0/nsample)
varia(:,:,0) = varia(:,:,0)/2.0
varia(:,:,nfrq-1) = varia(:,:,nfrq-1)/2.0
frq  = sdof@frq
splt = splt*conform(splt,frq,2)
vari1 = dim_sum_n(varia(:,0,:),1)

print(sdof)
printVarSummary(splt)
splt!2   = "frq"
splt&frq = log(frq)
splt@long_name = "power*frequency"
splt@units     = "1/day"

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.018
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
    
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "Period (day)"          ; xaxis
    resxy@tiYAxisString = "Power*frq"           ; yaxis
    ;resxy@tiYAxisString = "Variance"           ; yaxis
    ;resxy@trYLog  = True
    resxy@trYMinF = 0
    resxy@trYMaxF = max(splt(:,0,:))
    resxy@trXMinF = log(1.0/90.0) ;log(min(frq))
    resxy@trXMaxF = log(1.0/7.0)  ;log(max(frq)) 
    resxy@xyComputeXMin = True
    ;resxy@trXLog  = True
    resxy@tmXBMode    = "Explicit"
    resxy@tmXBValues  = log(1.0/p_draw)
    resxy@tmXBLabels  = p_draw
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case2 + " " + vari1
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@lgLabelFont            = font_index
    resxy@pmLegendWidthF   = 0.1
    resxy@pmLegendHeightF  = 0.25
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    if(ifrunave.ge.2) then 
    resxy@pmLegendOrthogonalPosF  = - 0.28   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = -0.28    ; move ref vector along X-axis
    else
    resxy@pmLegendOrthogonalPosF  = - 0.28   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = -0.28    ; move ref vector along X-axis
    end if
    
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

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    if(ifrunave.eq.0) then 
    resp@txString      = "sm="+sm+" pct="+pct+" no filter"
    end if
    if(ifrunave.eq.1) then 
    resp@txString      = "sm="+sm+" pct="+pct+" runave="+nave+"day"
    end if
    if(ifrunave.eq.2) then 
    resp@txString      = "sm="+sm+" pct="+pct+" bw_filter="+cb+"-"+ca+"day"
    end if
    if(ifrunave.eq.3) then 
    resp@txString      = "sm="+sm+" pct="+pct+" Fourier_filter="+cb+"-"+ca+"day"
    end if

if(allcase) then
    resxy@xyDashPatterns    = (/1,2,0,0,0,0/)
    resxy@xyLineThicknesses = (/3.0,3.0,3.0,5.0,5.0,5.0/)-1.5   ; make second line thicker
    resxy@xyLineColors      = (/"black","black","black","red","blue","grey"/)
    resxy@pmLegendDisplayMode = "Always"
    resxy@gsnLeftString   = pre_case(np) + ") " + month_s + "-" + month_e 
    plot(np) = gsn_csm_xy(wks, splt&frq, splt(:,0,:), resxy)     
    np = np + 1
;    gsn_panel(wks,plot,(/1,1/),resp)
else
    resxy@xyDashPatterns    = (/0,0,1,3/)
    resxy@xyLineThicknesses = (/2.0,3.0,3.0,3.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","red","red","red"/)
    resxy@pmLegendDisplayMode = "NoCreate"
    do nc = 0, ncase-1, 1
    resxy@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + vari1(nc)
    plot(np) = gsn_csm_xy(wks, splt&frq, splt(nc,:,:), resxy)     
    np = np + 1
    end do
    gsn_panel(wks,plot,(/1,ncase/),resp)
end if
end do
;gsn_panel(wks,plot,(/1,1/),resp)
;gsn_panel(wks,plot,(/2,2/),resp)
end
