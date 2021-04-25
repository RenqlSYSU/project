;***************************************************************************
;read preci, remove annual cycle, zonal average, 
;then yearly power spectrum is calculated
;make significance testing and draw the results
;
;                                           by Ql_Ren
;                                           2019/11/11
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql/project/2021variability"
area = (/"EA","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 10 , 22  , 30  /)
latn = (/ 50 , 30  , 45  /)
lonl = (/105 ,107  ,115  /)
lonr = (/125 ,130  ,130  /)
na = 0 ;"EA","TP","SEA","NEA"
lat_spcing = 5 
var_name = (/"preci","vort"/)
nv = 0 ;
lev = (/850,500,200,700,600,  100,10/)
nl = 0
stand = (/"","_stand"/)
ns = 0
month_s = 5 
month_e = 9
ifrunave= False;True ;
nave = 8 ;five days run ave

first_areaave = True ;False;
case_option  = 0 
diff = True;False;
perc = True;False
barshare = False;True ;
barvert  = False;True ;
landscap = False;True ;False mean portrait
siglvl = 0.05
p_draw = (/120,90,60,40,30,20,15,10,7,5/)

if(case_option.eq.0) then
    case1   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    if(diff) then 
    case2   = (/"CTRL","NG15S-CTRL","NG4550-CTRL","NGEN-CTRL"/) ;,"NG24h-CTRL"
    end if
    years   = ispan(1979,2005,1)   ;the calculate time
    path2   = "/home/ys17-19/renql"
    filein  = path2 +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                       "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                       "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                       "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
end if

fig_name = path + "/fig/" + case1(0) + "-" + month_s + month_e + stand(ns) + "_spec_" + var_name(nv)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

f    = addfile(filein(0)+years(9)+".daily.PRECC.nc","r")
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

spc  = new((/ncase,nyear,nfrq,nlat/),typeof(vars))
ave  = new((/ncase,indx ,nfrq,nlat/),typeof(vars))
pero = new((/ncase,indx-1,nfrq,nlat/),float) ;main period that pass the significance
prob = new((/ncase-1,nfrq,nlat/),float)
plot = new((/ncase/),graphic)
sig1 = new((/ncase/),graphic)
sig2 = new((/ncase/),graphic)
sig3 = new((/ncase/),graphic)

;---------------------------------------------------------------
;read the data, then calculate the yearly power spectra 
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 
  ;If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

do nc = 0, ncase-1 ,1
print("handle with " + case1(nc) + " " + area(na))

;========= read data ======================================
term := new((/nyear,ntime,nlat,nlon/),float)
if (nv.eq.0) then ;preci
do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    term(nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    term(nt,:,:,:) = term(nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
term = term*1000*3600*24 ;convert m/s to mm/day
end if

if (nv.eq.1) then ;vort
do nt = 0, nyear-1 ,1
    f  = addfile(filein(nc)+years(nt)+".daily.U.nc","r")
    u  = f->U(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    f  = addfile(filein(nc)+years(nt)+".daily.V.nc","r")
    v  = f->V(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    term(nt,:,:,:) = uv2vr_cfd(u,v,v&lat,v&lon,2) 
end do
term = term*10^5 ;10^-5 s^-1
end if

;========= data processing and calc spectrum ======================================
term  = linmsg_n(term, -1, 3)
term  = linmsg_n(term, -1, 2)
printVarSummary(term)
if (any(ismissing(term))) then
    print("term has missing values")
end if
opt  = 1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
dimh = 1
if(ifrunave) then 
term = runave_n_Wrap(term, nave, opt, dimh)
end if
term = term - conform(term,dim_avg_n(term,0),(/1,2,3/))  ;remove annual cycle
printVarSummary(term)

if(first_areaave) then 
    term := dim_avg_n_Wrap(term,3) ;zonal average
    if(ns.eq.1) then
    term = dim_standardize_n(term,1,1)
    end if
   
    do ny = 0, nlat-1, 1
    r1zsum = 0.0
    do nt = 0, nyear-1, 1
        ; calculate spectrum
        sdof := specx_anal(term(nt,:,ny),d,sm,pct) ;returns the degrees of freedom as a scalar.
        r1     = sdof@xlag1                        ; extract segment lag-1
        r1zsum = r1zsum  + 0.5*log((1+r1)/(1-r1)) ; sum the Fischer Z
        spc(nc,nt,:,ny) = sdof@spcx                ; sum spc of each segment
    end do 
    r1z  = r1zsum/nyear                 ; average r1z
    r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1); transform back, this is the mean r1
    spcavg  = dim_avg_n(spc(nc,:,:,ny),0)              ; average spectrum
    df      = 2.0*nyear   ; deg of freedom
    df@spcx = spcavg      ; assign the mean spc
    df@frq  = sdof@frq
    df@xlag1= r1          ; assign mean lag-1
    ave(nc,:,:,ny) = specx_ci(df, siglvl, 1.0-siglvl) 
    end do 

else

    spc1 = new((/nyear,nfrq,nlat,nlon/),typeof(term))
    ave1 = new((/indx ,nfrq,nlat,nlon/),typeof(term))
    if(ns.eq.1) then
    term = dim_standardize_n(term,1,1) ;dim_standardize_n(x, opt, dims)
    end if

    do ny = 0, nlat-1, 1
    do nx = 0, nlon-1, 1
    r1zsum = 0.0
    do nt = 0, nyear-1, 1
        ; calculate spectrum
        sdof := specx_anal(term(nt,:,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
        r1     = sdof@xlag1                        ; extract segment lag-1
        r1zsum = r1zsum  + 0.5*log((1+r1)/(1-r1)) ; sum the Fischer Z
        spc1(nt,:,ny,nx) = sdof@spcx                ; sum spc of each segment
    end do 
    r1z  = r1zsum/nyear                 ; average r1z
    r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1); transform back, this is the mean r1
    spcavg  = dim_avg_n(spc1(:,:,ny,nx),0)              ; average spectrum
    df      = 2.0*nyear   ; deg of freedom
    df@spcx = spcavg      ; assign the mean spc
    df@frq  = sdof@frq
    df@xlag1= r1          ; assign mean lag-1
    ave1(:,:,ny,nx) = specx_ci(df, siglvl, 1.0-siglvl) 
    end do 
    end do 
    ave(nc,:,:,:) = dim_avg_n(ave1,3) 
    spc(nc,:,:,:) = dim_avg_n(spc1,3) 
    term := dim_avg_n(term,3) 
end if
end do ;end case

frq = sdof@frq
ave = ave*conform(ave,frq,2)

print(sdof)
copy_VarMeta(vars(:,0),ave(0,0,0,:))
ave!2   = "frq"
ave&frq = log(frq)
ave@long_name = "power*frequency"
;ave@long_name = "frequency"
;ave@units     = "1/day"
printVarSummary(ave)

;---------------------------------------------------------------
;calc significance and difference, then draw the picture
;--------------------------------------------------------------
vari = dim_variance_n_Wrap(spc,1)
do nf = 0 ,ncase-2,1
prob(nf,:,:) = ttest(ave(0,0,:,:),vari(0,:,:),nyear,ave(nf+1,0,:,:),vari(nf+1,:,:),nyear,True,False)
end do
prob  = mask(prob,prob.lt.siglvl,True) ;to determine the dominant period
term1 = conform(pero,ave(:,0,:,:),(/0,2,3/))
pero  = mask(term1,term1.gt.ave(:,1:3,:,:),True)
delete([/term,term1,vari/])

if(diff) then 
do nf = 0 ,ncase-2,1
    if(perc) then 
    ave(nf+1,:,:,:) = ave(nf+1,:,:,:)/ave(0,:,:,:)*100.0-100.0
    else
    ave(nf+1,:,:,:) = ave(nf+1,:,:,:)-ave(0,:,:,:)
    end if
end do
end if

copy_VarMeta(ave(0,0,:,:),prob(0,:,:))
copy_VarMeta(ave(0,0,:,:),pero(0,0,:,:))

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 2.5
   font_index     = 22
   

wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    res@vpHeightF        = 0.8
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = True
    res@tiYAxisOn       = True
    res@tiXAxisFont     = font_index 
    res@tiYAxisFont     = font_index 
    res@tiYAxisString   = "day"
    res@tiXAxisString   = "Latitude (N)"

    res@tmYLMode        = "Explicit"
    res@tmYLValues      = log(1.0/p_draw)
    res@tmYLLabels      = p_draw

    res@tmXBMode        = "Manual"
    res@tmXBTickStartF  = lats(na)
    res@tmXBTickEndF    = latn(na)
    res@tmXBTickSpacingF= lat_spcing 

    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmYMinorGrid                   = True
    res@tmYMinorGridLineDashPattern    = 2
    res@tmYMinorGridLineColor          = "grey0"
    res@tmYMinorGridThicknessF         = 2.0
    
    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "grey0"
    res@tmYMajorGridThicknessF         = 2.0
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMajorOutwardLengthF = -0.
    res@tmYLMajorOutwardLengthF = -0.
    res@tmYRMajorOutwardLengthF = -0.
    
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    res@tmXBMinorOutwardLengthF = -0.
    res@tmYLMinorOutwardLengthF = -0.
    res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font

    res_sig                  = True       ; use plot options
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@gsnFrame         = False
    res_sig@gsnDraw          = False
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnLineThicknessF = 2
    ;res_sig@cnLineColor     = "black"
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.25
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "
    
    res_sig1 = res_sig
    res_sig1@cnFillPattern    = 7 
    res_sig2 = res_sig
    res_sig2@cnFillPattern    = 2 

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

if(nv.eq.0) then 
if(ns.eq.0) then 
    minlvl := (/30 ,30 ,30 ,30 /) 
    maxlvl := (/480,480,480,480/)
    spacig := (/30 ,30 ,30 ,30 /)
if(diff) then                  
    minlvl := (/30 ,-80,-80,-80/) 
    maxlvl := (/480, 80, 80, 80/) ;used for CTRL-NUDG
    spacig := (/30 , 10, 10, 10/)
end if                         
end if                         
                               
if(ns.eq.1) then               
    minlvl := (/0  ,0  ,0  ,0  /) 
    maxlvl := (/6.0,6.0,6.0,6.0/)
    spacig := (/0.4,0.4,0.4,0.4/)
if(diff) then                  
if(perc) then                  
    minlvl := (/0  ,-80,-80,-80/) 
    maxlvl := (/6.0, 80, 80, 80/) ;used for CTRL-NUDG
    spacig := (/0.4, 10, 10, 10/)
else
    minlvl := (/0  ,-0.8,-0.8,-0.8/) 
    maxlvl := (/6.0, 0.8, 0.8, 0.8/) ;used for CTRL-NUDG
    spacig := (/0.4, 0.1, 0.1, 0.1/)
end if
end if
end if
end if

if(nv.eq.1) then  ;vort
    minlvl := (/30 ,30 ,30 ,30 /) 
    maxlvl := (/480,480,480,480/)
    spacig := (/30 ,30 ,30 ,30 /)
if(diff) then              
    minlvl := (/30 ,-80,-80,-80/) 
    maxlvl := (/480, 80, 80, 80/) ;used for CTRL-NUDG
    spacig := (/30 , 10, 10, 10/)
end if
end if

np = 0
do nc = 0, ncase-1,1
    res@cnFillPalette := "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    if((nc.gt.0).and.diff) then
    res@cnFillPalette := "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
    end if 
    res@gsnLeftString   = pre_case(np) + ") " + case2(nc)
    if(lonl(na).eq.lonr(na)) then 
    res@gsnRightString  = lonl(na)+"N " + var_name(nv); +" ("+   
    else
    res@gsnRightString  = lonl(na)+"-"+lonr(na) +"E " + var_name(nv) 
    end if
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF      := spacig(nc)        ; set contour interval
    res@cnMaxLevelValF       := maxlvl(nc) ;32    ; set max contour level
    res@cnMinLevelValF       := minlvl(nc)       ; set min contour level
    plot(np)  = gsn_csm_time_lat(wks,ave(nc,0,:,:),res)
    print(max(ave(nc,0,:,:)))
    print(min(ave(nc,0,:,:)))

    if(diff.and.nc.gt.0) then 
        sig1(np) = gsn_csm_contour(wks,prob(nc-1,:,:),res_sig)
        overlay(plot(np),sig1(np))
    else
        ;sig1(np) = gsn_csm_contour(wks,pero(nc,1,:,:),res_sig1)
        ;overlay(plot(np),sig1(np))
        ;sig2(np) = gsn_csm_contour(wks,pero(nc,0,:,:),res_sig2)
        ;overlay(plot(np),sig2(np))
        sig3(np) = gsn_csm_contour(wks,pero(nc,0,:,:),res_sig)
        overlay(plot(np),sig3(np))
    end if
    np = np + 1
end do 

    resp@txString      = "sm="+sm+" pct="+pct+" runave="+ifrunave+" nave="+nave
    gsn_panel(wks,plot,(/1,ncase/),resp)
end

