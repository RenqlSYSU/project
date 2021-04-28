;***************************************************************************
;read var(nyear,nday),remove annual cycle and standardization
;then calc the wavelet for every year
;then average power,phase,signif,gws over all years.
;                                            by Ql_Ren
;                                           2018/03/30
;******************************************************************************
begin
path    = "/home/ys17-19/renql/project/2021variability"
area = (/"EA","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 27 , 22  , 30  /)
latn = (/ 37 , 30  , 45  /)
lonl = (/105 ,107  ,115  /)
lonr = (/125 ,130  ,130  /)
na = 0 ;"EA","TP","SEA","NEA"
var_name = (/"preci","vort"/)
nv = 0 ;
lev = (/850,500,200,700,600,  100,10/)
nl = 2
stand = (/"","_stand"/)
ns = 0
month_s = 1 
month_e = 12
ifrunave= False;True ;
nave = 8 ;five days run ave

first_areaave = True ;False;
case_option  = 0 
diff = False;True ;
perc = True ;False
barshare = True ;False;
barvert  = True ;False;
landscap = False;True ;False mean portrait
p_draw = (/120,90,60,40,30,20,15,10,7,5/)

if(case_option.eq.0) then
    case1   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    if(diff) then 
    case2   = (/"CTRL","NG15S-CTRL","NG4550-CTRL","NGEN-CTRL"/) ;,"NG24h-CTRL"
    end if
    years   = ispan(1979,1981,1)   ;the calculate time
    path2   = "/home/ys17-19/renql"
    filein  = path2 +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                       "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                       "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                       "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
end if

fig_name = path + "/fig/" + case1(0) + "-" + month_s + month_e + stand(ns) + "_wavelet_" + var_name(nv)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

f    = addfile(filein(0)+years(0)+".daily.PRECC.nc","r")
vars:= f->PRECC(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
ntime = dimsizes(select_time)
ncase = dimsizes(case1)
nyear = dimsizes(years)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)

;--the parameter for wavelet----------------
mother  = 0    ;0 = Morlet 1 = Paul 2 = DOG (derivative of Gaussian)
param   = 6.0  ;The mother wavelet parameter. If param < 0, then the default is used,
;Morlet k0 (wavenumber), default is 6. Paul m (order), default is 4. DOG m (m-th derivative), default is 2.

N       = ntime
dt      = 1    ;timesteps of time series, the unit here is day
s0      = dt   ;The smallest scale of the wavelet, which is typically equal to 2*dt. s0 = dt for Morlet; s0 = dt/4 for Paul
dj      = 0.25 ;The spacing between discrete scales, which is typically equal to 0.25. 
;A smaller value will give better scale resolution, but will be slower.
jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.))  
;The integer number of scales. Scales range from s0 up to s0*2^[(jtot-1)*dj].

npad    = N    ;The total number of points (including padding) to use for the wavelet transform
nadof   = 0    ;Currently ignored (set to zero).

noise   = 1    ;0 means use a white noise for significance test. 1 means use a red noise for significance test. (Most commonly, noise = 1.)
siglvl  = .05  ;The significance level to use
isigtest= 0    ; 0 means do a regular chi-square test, i.e. Eqn (18) from Torrence and Compo. 
;1 means do a "time-average" test on the global wavelet spectrum.

var  = new((/nyear,ntime/),float)
power  = new((/ncase,jtot,npad/),float)
sig    = new((/ncase,jtot,npad/),float)
gws    = new((/ncase,jtot/),float)
power1 = new((/jtot,npad/),float)
sig1   = new((/jtot,npad/),float)
gws1   = new((/jtot/),float)
plot   = new((/ncase/),graphic)
plot2  = new((/ncase/),graphic)

do nc = 0, ncase-1,1
print("start calc " + case1(nc) + var_name(nv) + lev(nl))
;--------------------------------------------------------------
;read the data
;--------------------------------------------------------------
if(nv.eq.0) then 
do nt = 0, nyear-1 ,1
    f    = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    term = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f    = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    term = term + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    var(nt,:) = wgt_areaave_Wrap(term,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
end do
var = var*1000*3600*24 ;convert m/s to mm/day
end if

if (nv.eq.1) then ;vort
do nt = 0, nyear-1 ,1
    f  = addfile(filein(nc)+years(nt)+".daily.U.nc","r")
    u  = f->U(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    f  = addfile(filein(nc)+years(nt)+".daily.V.nc","r")
    v  = f->V(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    var(nt,:) = uv2vr_cfd(u,v,v&lat,v&lon,2) 
end do
var = var*10^5 ;10^-5 s^-1
end if

opt  = 1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
dimh = 1
if(ifrunave) then 
var = runave_n_Wrap(var, nave, opt, dimh)
end if
var = var - conform(var,dim_avg_n(var,0),1)  ;remove annual cycle
if(ns.eq.1) then
var = dim_standardize_n(var,1,1) ;dim_standardize_n(x, opt, dims)
end if

power1 = 0
sig1   = 0
gws1   = 0
do nt = 0, nyear-1 ,1
    print("year " + years(nt))
    wave = wavelet(var(nt,:),mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
    power1 = power1 + onedtond(wave@power,(/jtot,npad/))
    sig1 = sig1 + power1/conform(power1,wave@signif,0)
    gws1 = gws1 + wave@gws  
    ;mean   = mean   + wave@mean
    ;stdev  = stdev  + wave@stdev
end do
power(nc,:,:) = power1/nyear
power(nc,:,:) = power(nc,:,:)/conform(power(nc,:,:),wave@scale ,0)   ;no bias, referenced by Liu et.al(2007)
sig(nc,:,:) = sig1/nyear
gws(nc,:) = gws1/nyear
end do;for case

power!1         = "period"
power&period    = wave@period
power!2         = "time"
power&time      = ispan(1,ntime,1) 
copy_VarMeta(power,sig)
sig = mask(sig,sig.ge.1.0,True)
printVarSummary(power)
printVarSummary(sig)

;==============================================================
;draw figure 
;==================================================================
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
    res@tiXAxisOn       = False
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
    res@tmYLValues          = p_draw ;ylabel 
    res@tmYLLabels          = p_draw ;ylabel
    res@tmLabelAutoStride   = True
    res@gsnYAxisIrregular2Log = True

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
    res@gsnRightString   = var_name(nv) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
    res@gsnCenterString  = "" 

    resp  = True    
    resp@gsnMaximize   = True
    if(landscap) then
    resp@gsnPaperOrientation = "landscape";"portrait"
    else
    resp@gsnPaperOrientation = "portrait"
    end if

    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn      = False   ; the label bar
    resp@gsnPanelLabelBar = True
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.7
    else
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.8 ;0.05
    resp@pmLabelBarHeightF = 0.05;0.8
    end if
    else
    resp@gsnPanelLabelBar = False; True
    res@lbLabelBarOn         = True   ; the label bar
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarZone               = 0
    if(barvert) then
    res@lbOrientation        = "Vertical";
    res@pmLabelBarWidthF             = 0.03
    res@pmLabelBarHeightF            = 0.17
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF          = 0.7
    res@pmLabelBarHeightF         = 0.08
    res@pmLabelBarOrthogonalPosF  = 0.05 ;y direction
    res@pmLabelBarParallelPosF    = 0.5
    end if
    end if
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if (nv.eq.0) then 
    minlvl := (/ 0 ,0 ,0 ,0 /) 
    maxlvl := (/ 7.5,7.5,7.5,7.5/) ;used for CTRL-NUDG
    spacig := (/ 0.5 ,0.5 ,0.5 ,0.5 /)
end if

do nc = 0,ncase-1,1
    res@cnFillPalette := "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    if((nc.gt.0).and.diff) then
    res@cnFillPalette := "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
    end if 
    res@gsnLeftString   = pre_case(nc) + ") " + case2(nc); + " " + vari1(nc)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF      := spacig(nc)        ; set contour interval
    res@cnMaxLevelValF       := maxlvl(nc) ;32    ; set max contour level
    res@cnMinLevelValF       := minlvl(nc)       ; set min contour level
    plot(nc)  = gsn_csm_contour(wks,power(nc,{2:120},:),res)
    plot2(nc) = gsn_csm_contour(wks,sig(nc,{2:120},:)  ,res_sig)
    overlay(plot(nc),plot2(nc))
    plot(nc) = ShadeCOI(wks, plot(nc), wave, power&time, False)
    print(max(power(nc,:,:)))
    print(min(power(nc,:,:)))
end do

gsn_panel(wks,plot,(/ncase,1/),resp)

end
