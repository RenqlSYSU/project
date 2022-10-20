;***************************************************************************
;read the spcx  
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/11/11
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql/project/TP_NUDG/spectral-20190120"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
diff = True
var_name = (/"U","V","T","vort"/)
area = (/"EA","TP","SEA","NEA","SNP","CNP","NNP"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 25 ,25  , 22  , 30  , 25  , 40  , 55  /)
latn = (/ 40 ,40  , 30  , 45  , 45  , 60  , 75  /)
lonl = (/ 0  ,70  ,107  ,115  ,140  ,140  ,140  /)
lonr = (/360 ,105 ,130  ,130  ,180  ,180  ,180  /)
nv = 3 ;"U","V","T","vort"
na = 4
month_s = 4
month_e = 9

if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h","NUDG24h"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG6h","NG24h"/) ;,"NG24h"
    if(diff) then 
    case2   = (/"CTRL","NG6h","NG6h-CTRL"/) ;,"NG24h"
    end if
    years   = ispan(1979,2005,1)   ;the calculate time
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case1   = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    case2    = (/"CTRL","noTP"/)
    years   = ispan(21,30,1)   ;the calculate time
end if

if(case_option.eq.2) then
    case1    = (/"F2000","NG6h"/)
    case2    = (/"F2000","NG6h"/)
    years    = ispan(1,30,1)   ;the calculate time
end if

lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
filein    = path + "/mdata/" + case1    + "-" + month_s + month_e + "_spec_" + var_name(nv) + "-" + area(na) + ".nc"
fig_name  = path + "/fig/"   + case1(0) + "-" + month_s + month_e + "_spec_" + var_name(nv) + "_lev_frq"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

f    = addfile(filein(0),"r")
spcx = f->spcx
frq  = f->frq
nfrq = dimsizes(spcx&frq)
nlev = dimsizes(spcx&lev)
ncase = dimsizes(case1)
var  = new((/ncase,nfrq ,nlev/),float)
plot = new((/ncase*3/),graphic)

np = 0
wks = gsn_open_wks(fig_out,fig_name)  

do na = 4,6,1
filein    = path + "/mdata/" + case1 + "-" + month_s + month_e + "_spec_" + var_name(nv) + "-" + area(na) + ".nc"
do nc = 0, ncase-1, 1
    f    = addfile(filein(nc),"r")
    var(nc,:,:) = f->spcx
end do
printVarSummary(var)
copy_VarMeta(spcx,var(0,:,:))

if(nv.eq.3) then 
    var = var*10^11
end if

if(diff) then 
    var(2,:,:) = var(1,:,:)/var(0,:,:)-1.0
end if

if(nfrq.gt.200) then
print(frq(0:1000:100))
print(var(0,0:1000:100,5))
nfrq2 = 100
nspc  = toint(nfrq/nfrq2)
nsta  = nfrq - nspc*nfrq2
print("nfrq = " + nfrq2 + ", spacing = " + nspc + ", start point = " + nsta)
var2 = new((/ncase,nfrq2,nlev/),float)
frq2 = new((/nfrq2/),float)
calc_option = 1
if(calc_option.eq.0) then
    opt = 0 ;Reserved for future use. It is currently not used, but set it to 0
    cyclic = False ;indicate whether the dim-th dimension of data is cyclic
    frq2 = fspan(0.001,0.5,nfrq2)
    var2 = linint1_n_Wrap(frq,var,cyclic,frq2,opt,1)
else
    do nfr = 0, nfrq2-1, 1
        n1 = nsta - 1 + nspc * nfr
        n2 = nsta - 1 + nspc *(nfr+1)
        n3 = (n1+n2)/2
        var2(:,nfr,:) = dim_avg_n_Wrap(var(:,n1:n2,:),1)
        frq2(nfr) = frq(n3)
    end do
    var2!1   = "frq"
    var2&frq = frq2
end if
else
var2 = var
frq2 = frq
end if
var2!0 = "case"
printVarSummary(var2)
print(frq2(0:20))
print(var2(0,0:20,5))

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
   mp_thick       = 5.0
   font_index     = 22
   
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    res@vpHeightF        = 0.6
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = True
    res@tiYAxisOn       = True
    res@tiXAxisFont     = font_index 
    res@tiYAxisFont     = font_index 
    res@tiXAxisString   = "Frequency"

    res@tmXBMode        = "Explicit"
    res@tmXBValues      = fspan(0.0,0.5,6)
    res@tmXBLabels      = fspan(0.0,0.5,6)
    
    res@gsnPresHgtHeightLabelOn = False
    res@tiYAxisOn = False
    res@tmYRLabelsOn = False

    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmYMinorGrid                   = False
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
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Vertical";"Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2
    res@pmLabelBarWidthF  = 0.1
    res@pmLabelBarHeightF = 0.6

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString  = var_name(nv) + " " + area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = False
    resp@lbOrientation    = "Horizontal"
    resp@lbLabelStride    = 2
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.05
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    minlvl := (/10 ,10 ,10 /) ;used for East Asian moist transiation
    maxlvl := (/170,170,170/) ;used for CTRL-NUDG
    spacig := (/10 ,10 ,10 /)
if(diff) then 
    minlvl := (/10 ,10 ,-0.8/) ;used for East Asian moist transiation
    maxlvl := (/490,490, 0.8/) ;used for CTRL-NUDG
    spacig := (/30 ,30 , 0.1/)
end if

;ip   = ind(frq.ge.0.001)    ; all indices for "short" periods

do nc = 0, ncase-1,1
    res@cnFillPalette := "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
if((nc.eq.2).and.diff) then
    res@cnFillPalette := "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
end if 
    res@gsnLeftString   = pre_case(np) + ") " + case2(nc)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF      := spacig(nc)        ; set contour interval
    res@cnMaxLevelValF       := maxlvl(nc) ;32    ; set max contour level
    res@cnMinLevelValF       := minlvl(nc)       ; set min contour level
    plot(np)  = gsn_csm_pres_hgt(wks,var2(case|nc,lev|:,frq|:),res)
    print(max(var2(nc,:,:)))
    print(min(var2(nc,:,:)))
    np = np + 1
end do
delete(var2)
end do 

    gsn_panel(wks,plot,(/3,ncase/),resp)
end

