;***************************************************************************
;read the spcx  
;then do the spectral analysis and plot the figure
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
lat_spcing = 5 
var_name = (/"preci"/)
stand = (/"","_stand"/)
ns = 1
month_s = 5 
month_e = 8

case_option  = 0 ;0 means NUDG experiment; 1 means noTP experiment
diff = True;False;
perc = True;False
barshare = False;True ;
barvert  = False;True ;
landscap = False;True ;False mean portrait
na = 0 ;"EA","TP","SEA","NEA"
nv = 0 ;"U","V","T","vort"
siglvl = 0.05
p_draw = (/120,90,60,40,30,20,15,10,7,5/)

if(case_option.eq.0) then
    case1   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    if(diff) then 
    case2   = (/"CTRL","NG15S-CTRL","NG4550-CTRL","NGEN-CTRL"/) ;,"NG24h-CTRL"
    end if
    years   = ispan(1979,2005,1)   ;the calculate time
end if

filein   = path + "/data/" + case1 + "-" + month_s + month_e + stand(ns) + "_spec_" + var_name(nv) + "-" + area(na) + ".nc"
fig_name = path + "/fig/" + case1(0) + "-" + month_s + month_e + stand(ns) + "_spec_" + var_name(nv)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

f    = addfile(filein(0),"r")
spcx = f->spcx(:,:,:,{lats(na):latn(na)},{lonl(na):lonr(na)})
frq  = f->frq
nfrq = dimsizes(spcx&frq)
nlon = dimsizes(spcx&lon)
nlat = dimsizes(spcx&lat)
ncase = dimsizes(case2)
nyear = dimsizes(years)
indx = 4 ;the result of specx_ci
var  = new((/ncase,nyear,indx  ,nfrq,nlat,nlon/),float)
pero = new((/ncase,indx-1,nfrq,nlat/),float)
prob = new((/ncase-1,indx,nfrq,nlat/),float)
plot = new((/ncase/),graphic)
sig1 = new((/ncase/),graphic)
sig2 = new((/ncase/),graphic)
sig3 = new((/ncase/),graphic)

do nc = 0, ncase-1, 1
    f    = addfile(filein(nc),"r")
    var(nc,:,:,:,:,:) = f->spcx(:,:,:,{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
printVarSummary(var)
var := dim_avg_n_Wrap(var,5)
print(var(0, 0, 0, 0:50:10, 5))

vari = dim_variance_n_Wrap(var,1)
var := dim_avg_n_Wrap(var,1)
do nf = 0 ,ncase-2,1
prob(nf,:,:,:) = ttest(var(0,:,:,:),vari(0,:,:,:),nyear,var(nf+1,:,:,:),vari(nf+1,:,:,:),nyear,True,False)
end do
term = conform(pero,var(:,0,:,:),(/0,2,3/))
pero = mask(term,term.gt.var(:,1:3,:,:),True)
prob = mask(prob,prob.lt.siglvl,True)
delete([/term,vari/])

if(diff) then 
do nf = 0 ,ncase-2,1
    if(perc) then 
    var(nf+1,:,:,:) = var(nf+1,:,:,:)/var(0,:,:,:)*100.0-100.0
    else
    var(nf+1,:,:,:) = var(nf+1,:,:,:)-var(0,:,:,:)
    end if
end do
end if

if(nfrq.gt.200) then 
nfrq2 = 100
nspc  = toint(nfrq/nfrq2)
nsta  = nfrq - nspc*nfrq2
print("nfrq = " + nfrq2 + ", spacing = " + nspc + ", start point = " + nsta)
var2 = new((/ncase,nfrq2,nlev,nlon/),float)
frq2 = new((/nfrq2/),float)
calc_option = 1
if(calc_option.eq.0) then 
    frq2 = fspan(0.001,0.5,nfrq2)
    var2 = linint1_n_Wrap(frq,var,False,frq2,0,1)
else
    do nfr = 0, nfrq2-1, 1
        n1 = nsta - 1 + nspc * nfr
        n2 = nsta - 1 + nspc *(nfr+1)
        n3 = (n1+n2)/2
        var2(:,nfr,:,:) = dim_avg_n_Wrap(var(:,n1:n2,:,:),1)
        frq2(nfr) = frq(n3)
    end do
    var2!1   = "frq"
    var2&frq = frq2
end if
printVarSummary(var2)
print(frq2(0:20))
print(frq(0:100:10))
print(var2(0,0:20,0,5))
print(var(0,0:100:10,0,5))
else
var2 = var
frq2 = frq
end if

var2&frq = 1.0/frq2
copy_VarMeta(var2(0,0,:,:),prob(0,0,:,:))
copy_VarMeta(var2(0,0,:,:),pero(0,0,:,:))

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
    res@tiYAxisString   = "Frequency (1/day)"
    res@tiXAxisString   = "Latitude (N)"

    res@tmYLMode        = "Explicit"
    res@tmYLValues      = p_draw
    res@tmYLLabels      = p_draw
    ;res@tmYLValues      = fspan(0.0,0.5,6)
    ;res@tmYLLabels      = fspan(0.0,0.5,6)
    res@gsnYAxisIrregular2Log = True

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
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/10,20,40,80,120,  160,200,300,400,600,  800,1000,1200,1400,1600,  1800,2000/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 10         ; set min contour level
    ;res@cnMaxLevelValF        = 330       ; set max contour level
    ;res@cnLevelSpacingF       = 20        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
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
    plot(np)  = gsn_csm_time_lat(wks,var2(nc,0,:,:),res)
    print(max(var2(nc,0,:,:)))
    print(min(var2(nc,0,:,:)))

    if(diff.and.nc.gt.0) then 
        sig1(np) = gsn_csm_contour(wks,prob(nc-1,0,:,:),res_sig)
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

    gsn_panel(wks,plot,(/1,ncase/),resp)
end

