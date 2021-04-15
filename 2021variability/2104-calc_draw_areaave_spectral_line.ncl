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
area = (/"EA","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 27 , 22  , 30  /)
latn = (/ 37 , 30  , 45  /)
lonl = (/105 ,107  ,115  /)
lonr = (/125 ,130  ,130  /)
var_name = (/"preci"/)
stand = (/"","_stand"/)
ns = 0
month_s = 5 
month_e = 9

first_areaave = True ;False;
case_option  = 0 
siglvl = 0.1
na = 0 ;"EA","TP","SEA","NEA"
nv = 0 ;"U","V","T","vort"
p_draw = (/120,90,60,40,30,20,15,10,7,5/)

if(case_option.eq.0) then
    case1   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    years   = ispan(1979,2005,1)   ;the calculate time
    path2   = "/home/ys17-19/renql"
    filein  = path2 +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                       "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                       "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                       "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
end if

fig_name = path + "/fig/" + case1(0) + "-" + month_s + month_e + stand(ns) + "_area_spec_line_" + var_name(nv)
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
nfrq = toint(ntime/2.0) 
indx = 4 ;the result of specx_ci

var  = new((/ncase,nyear,indx ,nfrq/),float)
ave  = new((/ncase      ,indx ,nfrq/),float)
plot = new((/ncase/),graphic)
vari2 = new((/ncase/),float)
vari3 = new((/ncase/),float)

;---------------------------------------------------------------
;read the data, then calculate the yearly power spectra 
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0.1;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

do nc = 0, ncase-1 ,1
print("handle with " + case1(nc) + " " + area(na))
term := new((/nyear,ntime,nlat,nlon/),float)

do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    term(nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    term(nt,:,:,:) = term(nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
term  = linmsg_n(term, -1, 3)
term  = linmsg_n(term, -1, 2)
if (any(ismissing(term))) then
    print("term has missing values")
end if
term = term*1000*3600*24 ;convert m/s to mm/day
nave = 5 ;five days run ave
opt  = 1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
dimh = 1
term = runave_n_Wrap(term, nave, opt, dimh)
term = term - conform(term,dim_avg_n(term,0),(/1,2,3/))  ;remove annual cycle
printVarSummary(term)

if(first_areaave) then 
    term := wgt_areaave_Wrap(term,1.0,1.0,0) ;area average
    if(ns.eq.1) then
    term = dim_standardize_n(term,1,1) ;dim_standardize_n(x, opt, dims)
    end if

    do nt = 0, nyear-1, 1
    ; calculate spectrum
    sdof := specx_anal(term(nt,:),d,sm,pct) ;returns the degrees of freedom as a scalar.
    var(nc,nt,:,:) = specx_ci(sdof, siglvl, 1.0-siglvl) 
    end do 
    sdof := specx_anal(dim_avg_n(term,0),d,sm,pct) ;returns the degrees of freedom as a scalar.
    ave(nc,:,:) = specx_ci(sdof, siglvl, 1.0-siglvl) 

else

    var1  = new((/nyear,indx ,nfrq,nlat,nlon/),float)
    ave1  = new((/      indx ,nfrq,nlat,nlon/),float)
    if(ns.eq.1) then
    term = dim_standardize_n(term,1,1) ;dim_standardize_n(x, opt, dims)
    end if

    do ny = 0, nlat-1, 1
    do nx = 0, nlon-1, 1
    sdof := specx_anal(dim_avg_n(term(:,:,ny,nx),0),d,sm,pct) ;returns the degrees of freedom as a scalar.
    ave1(:,:,ny,nx) = specx_ci(sdof, siglvl, 1.0-siglvl) 
    do nt = 0, nyear-1, 1
    ; calculate spectrum
    sdof := specx_anal(term(nt,:,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
    var1(nt,:,:,ny,nx) = specx_ci(sdof, siglvl, 1.0-siglvl) 
    end do 
    end do 
    end do 
    var(nc,:,:,:) = wgt_areaave_Wrap(var1,1.0,1.0,0) ;area average
    ave(nc,:,:)   = wgt_areaave_Wrap(ave1,1.0,1.0,0) ;area average
end if
vari2(nc) = dim_avg_n(dim_variance_n(term,1),0)
vari3(nc) = dim_variance(ndtooned(term))
end do ;end case

var := dim_avg_n_Wrap(var,1)
ave = var
;ave(:,0,:) = var(:,0,:)
ave = ave*(1.0/ntime)
ave(:,:,0) = ave(:,:,0)/2.0
ave(:,:,nfrq-1) = ave(:,:,nfrq-1)/2.0
frq = sdof@frq
;ave = ave*conform(ave,frq,2)

;calc total variance
vari1 = dim_sum_n(ave(:,0,:),1)
print(vari1)
print(vari2)
print(vari3)

print(sdof)
printVarSummary(ave)
ave!2   = "frq"
ave&frq = frq
ave@long_name = "frequency"
ave@units     = "1/day"

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

wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file

    resxy     = True
    ;resxy@vpWidthF  = 0.45
    ;resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = (/0,0,0,0/)
    resxy@xyLineThicknesses = (/5.0,3.0,3.0,3.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"grey","red","blue","black"/)
    ;resxy@xyDashPatterns    = (/0,0,1,3/)
    ;resxy@xyLineThicknesses = (/2.0,3.0,3.0,3.0/)   ; make second line thicker
    ;resxy@xyLineColors      = (/"black","red","red","red"/)
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "Period (day)"          ; xaxis
    resxy@tiYAxisString = "Variance"           ; yaxis
    ;resxy@trYLog  = True
    ;resxy@trYMinF = 0
    ;resxy@trYMaxF = 4.0 
    resxy@xyComputeXMin = True
    resxy@trXLog  = True
    resxy@tmXBMode    = "Explicit"
    resxy@tmXBValues  = 1.0/p_draw
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
    resxy@pmLegendHeightF  = 0.18
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = - 0.35   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.25    ; move ref vector along X-axis
    
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
    resxy@gsnRightString   = var_name(nv) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"

;ip   = ind(frq.le.0.2)    ; all indices for "short" periods
ip   = ind(frq.ge.0.01)    ; all indices for "short" periods
np = 0
;do nc = 0, ncase-1, 1
    resxy@gsnLeftString   = pre_case(np) + ") " ;+ case2(nc)
    plot(np) = gsn_csm_xy(wks, ave&frq, ave(:,0,:), resxy)     
    np = np + 1
;end do

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

;    resp@txString      = fig_title
    ;gsn_panel(wks,plot,(/1,ncase/),resp)
    gsn_panel(wks,plot,(/1,1/),resp)
end