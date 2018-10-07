******************************************************************************
;filename:      draw_1p4x3_month_pc_wavelet.ncl
;input_file:    var_name(nv) + area + _EOFn3.nc
;----------------------------------------------------------------------
; Concepts illustrated:
;   - read EOF time series data and eval
;   - calc wavelet of EOF time series
;   - Draw the
;                                            by Ql_Ren
;                                           2018/10/07
;*********************************************************************************
begin
unit     = (/"~S~o~N~C","hPa","W/m2","W/m2"/)    ;this can be used as right string 
var_name = (/"TS","PSL","FLNT","FSNT"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)
nv = 0

;-----EOF area-------------- 
lats = -30
latn = 30
lonl = 120
lonr = 260
;lev  = 200
neval = 4
lat_spacing = 15
lon_spacing = 20

;----path of the file and output file-----------------------------
path      = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename  = path + "/mdata/" + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E_DJF_covEOF" + neval + ".nc"
fig_name  = path + "/fig/"   + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E_DJF_covEOF_wavelet" + neval
fig_out   = "pdf";"x11"

;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
ts_var  = dim_standardize_n_Wrap(f->eof_ts,1,2)    ;(ncase,neval,nyear) 
right_title   = f->pcvar   ;the percent variance associated with each eigenvalue.
center_title  = f->eval    ;eigenvalues
printVarSummary(ts_var)
printVarSummary(right_title)

ncase = dimsizes(ts_var&case)
ntime = dimsizes(ts_var&time)
gap_time = toint(ntime/30)
sig1  = new((/ncase,neval/),logical) 
plot  = new((/ncase*(neval-1)/),graphic)
plot2 = new((/ncase*(neval-1)/),graphic)

do nc = 0,ncase-1,1
sig1(nc,:) = eofunc_north(right_title(nc,:),ntime,True)
end do
subtitle_color = where(sig1,"red","black")

;---- compute wavelet----------------------------------
mother  = 0    ;0 = Morlet, 1 = Paull, 2 = DOG (derivative of Gaussian)
param   = 6.0  ;The mother wavelet parameter. If param < 0, then the default is used,
;Morlet k0 (wavenumber), default is 6. Paul m (order), default is 4. DOG m (m-th derivative), default is 2.

N       = ntime
dt      = 1    ;timesteps of time series, the unit here is month
s0      = dt   ;The smallest scale of the wavelet, which is typically is equal to 2*dt. s0 = dt for Morlet; s0 = dt/4 for Paul
dj      = 0.25 ;The spacing between discrete scales, which is typically equal to 0.25. A smaller value will give better scale resolution, but will be slower.
jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.))  
;The integer number of scales. Scales range from s0 up to s0*2^[(jtot-1)*dj].

npad    = N    ;The total number of points (including padding) to use for the wavelet transform
nadof   = 0    ;Currently ignored (set to zero).

noise   = 1    ;0 means use a white noise for significance test. 1 means use a red noise for significance test. (Most commonly, noise = 1.)
siglvl  = .05  ;The significance level to use
isigtest= 0    ; 0 means do a regular chi-square test, i.e. Eqn (18) from Torrence and Compo. 1 means do a "time-average" test on the global wavelet spectrum.

power = new((/ncase, neval, jtot, npad/),float)
sig2  = new((/ncase, neval, jtot, npad/),float)
do nc = 0, ncase-1, 1
do ne = 0, neval-1, 1
    wave = wavelet(ts_var(nc,ne,:),mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
    printVarSummary(wave)
    power(nc,ne,:,:) = onedtond (wave@power,(/jtot,npad/) )
    sig2(nc,ne,:,:)  = power(nc,ne,:,:)/conform(power(nc,ne,:,:),wave@signif,0) ;compute significance ( >= 1 is significant)
end do
end do

power!2         = "period"
power&period    = wave@period
power!3         = "time"
power&time      = ts_var&time
copy_VarMeta(power,sig2)
sig2 = mask(sig2,sig2.ge.1.0,True)
printVarSummary(power)
printVarSummary(sig2)

exponent = ispan(1,6,1)
ylabel   = exp( exponent*log(2.0) )
print(ylabel)

;==============================================================
;draw figure 
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.02
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@vpHeightF        = 0.15       ; height and width of plot
    res@vpWidthF         = 0.45

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
    res@tiXAxisOn       = False ;True
    res@tiXAxisString   = "Year"
    
    res@tiYAxisOn       = False
    res@tiYAxisFont     = font_index
    res@tiYAxisFontHeightF = subtitle_font 
    res@tiYAxisString   = "period(months)"
    
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = ispan(1,ntime,gap_time)
    res@tmXBLabels      = ispan(1,30,1)
    
    res@trYReverse          = True   ; reverse y-axis
    res@tmYLMode            = "Explicit"
    res@tmYLValues          = ylabel 
    res@tmYLLabels          = ylabel
    res@tmLabelAutoStride   = True

    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "black"
    res@tmYMajorGridThicknessF         = 4.0

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
    res@cnInfoLabelOn    = False
    res@cnFillMode          = "RasterFill"          ; turn on raster mode
    res@cnRasterSmoothingOn = True                  ; turn on raster smoothing

    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = cnlevels
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 0.5       ; set min contour level
    res@cnMaxLevelValF        = 7.5       ; set max contour level
    res@cnLevelSpacingF       = 0.5        ; set contour interval

;    res@gsnSpreadColors      = True              ; turn of the full color map
;    res@gsnSpreadColorStart  = 2                ; start at color 17
;    res@gsnSpreadColorEnd    = 21                ; end at color 200

    res@lbLabelBarOn         = False  ; the label bar
    res@lbOrientation        = "Horizontal"   ;"Vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-2
    
;----------Set res sig-----------
    res_sig                  = True
    res_sig@gsnDraw          = False
    res_sig@gsnFrame         = False
    
    res_sig@tiXAxisOn       = False
    res_sig@tiYAxisOn       = False
    res_sig@tmXBMode        = "Explicit"
    res_sig@tmXBValues      = ispan(1,ntime,gap_time)
    res_sig@tmXBLabels      = ispan(1,30,1)
    res_sig@trYReverse          = True   ; reverse y-axis
    res_sig@tmYLMode            = "Explicit"
    res_sig@tmYLValues          = ylabel
    res_sig@tmYLLabels          = ylabel
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
    res@gsnCenterString   = "" ;draw_var + " (" + unit +")"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Horizontal" ;"Vertical" ;"Horizontal"
    resp@lbLabelStride    = 1
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@lbLabelFontHeightF   = tick_font              ; make labels smaller
    resp@lbBoxLineThicknessF  = mp_thick-2
    resp@gsnPaperOrientation  = "landscape"

do nc = 0,ncase-1,1
do ne = 0,neval-2,1
    np = nc*(neval-1) + ne
    res@gsnLeftString    = case(nc)
    res@gsnCenterString  = "PC " + (ne+1)
    res@gsnRightStringFontColor    = subtitle_color(nc,ne)
    res@gsnRightString   = center_title(nc,ne)+"  "+right_title(nc,ne)+"%" 
    plot(np)  = gsn_csm_contour(wks,power(nc,ne,{1:64},:),res)
    plot2(np) = gsn_csm_contour(wks, sig2(nc,ne,{1:64},:),res_sig)
    overlay(plot(np),plot2(np))
    plot(np) = ShadeCOI(wks, plot(np), wave, ts_var&time, False)
    print(max(power(nc,ne,:,:)))
    print(min(power(nc,ne,:,:)))
end do
end do

gsn_panel(wks,plot,(/ncase,(neval-1)/),resp)
end

