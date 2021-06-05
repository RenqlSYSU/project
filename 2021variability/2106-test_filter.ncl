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
lats = (/ 27, 32 , 25 , 40 /)
latn = (/ 37, 32 , 25 , 40 /)
lonl = (/105, 118, 115, 115/)
lonr = (/125, 118, 115, 115/)
var_name = (/"preci","vort"/)
unit = (/"mm/day","s~N~-1~S~10~N~-5~S~"/)
nv = 0 ;
na = 0
lev = (/850,500,200,700,600,  100,10/)
nl = 1
stand = (/"","_stand"/)
ns = 0
nave = 8 ;five days run ave
ca = 80.0;90.0 ;day
cb = 8.0 ;2.5  ;10.0 ;
month_s = 6 
month_e = 8

if(nv.eq.0) then
    case1   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"CN05.1","APHRO","CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    years   = ispan(1985,1989,1)   ;the calculate time
    path2   = "/home/ys17-19"
    filein  = path2 +(/"/data-observation/CN05.1/CN05.1_Pre_1961_2018_daily_025x025.nc",\
                       "/data-observation/APHRO_MA_025deg_V1003R1/APHRO_MA_025deg_V1003R1.",\; + years(nt) + ".nc"
                       "/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                       "/renql/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
end if
fig_name = path + "/fig/" + case1(0) + "-" + month_s + month_e + "_filter_line_" + var_name(nv)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

f    = addfile(filein(3)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
ntime = dimsizes(select_time)
ncase = dimsizes(case1)
nyear = dimsizes(years)
narea = dimsizes(lats)

legend = (/"raw",nave+" runave",cb+"-"+ca+" bwfilter"/)
term   = new((/nyear,ntime/),float)
plot   = new((/ncase/),graphic)
vari2  = new((/ncase/),float)
vari3  = new((/ncase/),float)

;---------------------------------------------------------------
;read the data, then calculate the yearly power spectra 
;--------------------------------------------------------------
wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file
do na = 0,narea-1,1
var  := new((/ncase,3,nyear,ntime/),float)
do nc = 0, ncase-1 ,1
print("")
print("handle with " + case1(nc))

;========= read data ======================================
if(nv.eq.0) then ;preci
if(nc.eq.0) then
    f   = addfile(filein(nc),"r")
    t   = cd_calendar(f->time,0)
    select_time1 = ind(t(:,0).ge.years(0).and.t(:,0).le.years(nyear-1) \
                    .and.t(:,1).ge.month_s.and.t(:,1).le.month_e)
    term = reshape(wgt_areaave(f->pre(select_time1,{lats(na):latn(na)},{lonl(na):lonr(na)}),1.0,1.0,0),(/nyear,ntime/));mm/day
end if
if(nc.eq.1) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".nc","r")
    if(isleapyear(years(nt))) then ;leap year,366 day
    term(nt,:) = wgt_areaave(f->precip(select_time+1,{lats(na):latn(na)},{lonl(na):lonr(na)}),1.0,1.0,0)
    else
    term(nt,:) = wgt_areaave(f->precip(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}),1.0,1.0,0);mm/day
    end if
    end do
end if
if(nc.ge.2) then
    do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    term(nt,:) = wgt_areaave(f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}) ,1.0 ,1.0 ,0)
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    term(nt,:) = term(nt,:) + wgt_areaave(f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}) ,1.0 ,1.0 ,0)
    end do
    term = term*1000*3600*24 ;convert m/s to mm/day
end if
end if

;-------------------vort------------------------------
if(nv.eq.1) then ;vort
do nt = 0, nyear-1 ,1
    f  = addfile(filein(nc)+years(nt)+".daily.U.nc","r")
    u := f->U(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    f  = addfile(filein(nc)+years(nt)+".daily.V.nc","r")
    v := f->V(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    term(nt,:) = wgt_areaave(uv2vr_cfd(u,v,v&lat,v&lon,2)*10^5 ,1.0 ,1.0 ,0)
end do
end if

;========= data processing and calc spectrum ======================================
var(nc,0,:,:) = term

opt  = -1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
dimh = 1
var(nc,1,:,:) = runave_n_Wrap(term, nave, opt, dimh)

fca  = 1.0/ca
fcb  = 1.0/cb
opt1 = True
opt1@remove_mean = True
annu = conform(term,dim_avg_n(term,0),(/1/))  ;remove annual cycle
var(nc,2,:,:) = bw_bandpass_filter((term-annu),fca,fcb,opt1,1)+annu
;fcf = fspan(0, 0.5, ntime/2) 
;filt_cf = ind(fcf.lt.fca.or.fcf.gt.fcb)
;cf = ezfftf_n(term,1)
;cf(:,:,filt_cf) = 0.0  
;var(nc,2,:,:) = ezfftb_n(cf,cf@xbar,2)
end do
vari := dim_variance_n((var-conform(var,dim_avg_n(var,2),(/0,1,3/))),3)
vari := dim_avg_n(vari,2)
var := reshape(var,(/ncase,3,nyear*ntime/))
;vari = dim_variance_n(var,2)
;vari = dim_avg_n(var*var,2)
xtim = ispan(0,nyear*ntime-1,1)

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
    resxy@tiYAxisString = var_name(nv)+"("+unit(nv)+")"           ; yaxis
    resxy@tiXAxisFontHeightF = tick_font
    resxy@tiYAxisFontHeightF = tick_font
    resxy@tiXAxisFont = font_index
    resxy@tiYAxisFont = font_index
    resxy@trYMinF = min(var)
    resxy@trYMaxF = max(var)
    resxy@trXMinF = 0 
    resxy@trXMaxF = nyear*ntime-1
    resxy@xyComputeXMin = True
    resxy@tmXBMode    = "Explicit"
    resxy@tmXBValues  = ispan(0,nyear*ntime-1,ntime) 
    resxy@tmXBLabels  = years
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = legend 
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
    if(nv.eq.1) then 
    resxy@gsnRightString   = lev(nl) + "hPa " + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
    else
    resxy@gsnRightString   =  " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
    end if

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    resxy@xyDashPatterns    = (/0,0,0/)
    resxy@xyLineThicknesses = (/2.0,2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","red","blue"/)
    
    np = 0
    do nc = 0, ncase-1, 1
    resxy@gsnLeftString   = pre_case(np) + ") " + case2(nc)+": "+vari(nc,0)+"; "+vari(nc,1)+"; "+vari(nc,2)
    plot(np) = gsn_csm_xy(wks, xtim, var(nc,:,:), resxy)     
    np = np + 1
    end do
    gsn_panel(wks,plot,(/ncase,1/),resp)
end do
end
