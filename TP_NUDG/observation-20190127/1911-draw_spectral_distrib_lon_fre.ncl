;***************************************************************************
;read the spcx  
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/11/11
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127"
area = (/"EA","TP","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 35 ,25  , 22  , 30  /)
latn = (/ 35 ,40  , 30  , 45  /)
lonl = (/ 0  ,70  ,107  ,115  /)
lonr = (/360 ,105 ,130  ,130  /)
lon_spcing = 45
lev  = (/500,400,300,200/)
var_name = (/"U","V","T","vort"/)

area_option = True
diff = True
na = 0 ;"EA","TP","SEA","NEA"
nv = 3 ;"U","V","T","vort"
nl = 3 ;850,500,400,300,200
case = "ERA-Interim"

if(area_option) then
filein   = path + "/mdata/" + case + "_allyear_spec_" + var_name(nv) + "-" + area(na) + "-area.nc"
fig_name = path + "/fig/" + case + "_allyear_spec_" + var_name(nv) + "-" + lats(na) + "N_" + lev(nl)
else
filein   = path + "/mdata/" + case + "_allyear_spec_" + var_name(nv) + "-" + area(na) + ".nc"
fig_name = path + "/fig/" + case + "_allyear_spec_" + var_name(nv) + "-" + area(na)
end if
pre_case  = (/"(c","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_title = "Power Spectrum"
fig_out   ="pdf";"x11"

f    = addfile(filein,"r")
if(area_option) then 
spcx = f->spcx(:,{lev(nl)},{lats(na)},{lonl(na):lonr(na)})
else
spcx = f->spcx(:,{lev(nl)},{lonl(na):lonr(na)})
end if
frq  = f->frq
nfrq = dimsizes(spcx&frq)
nlon = dimsizes(spcx&lon)
var  = new((/nfrq ,nlon/),float)
plot = new((/6/),graphic)

   wks = gsn_open_wks(fig_out,fig_name)  
  setvalues NhlGetWorkspaceObjectId()
      "wsMaximumSize" : 300000000
  end setvalues
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue

do nl=0,3,1
f    = addfile(filein,"r")
if(area_option) then 
var = f->spcx(:,{lev(nl)},{lats(na)},{lonl(na):lonr(na)})
else
var = f->spcx(:,{lev(nl)},{lonl(na):lonr(na)})
end if
printVarSummary(var)
print(var(0:1000:100,5))

if(nv.eq.3) then 
    var = var*10^11
end if

nfrq2 = 100
nspc  = toint(nfrq/nfrq2)
nsta  = nfrq - nspc*nfrq2
print("nfrq = " + nfrq2 + ", spacing = " + nspc + ", start point = " + nsta)
var2 = new((/nfrq2,nlon/),float)
frq2 = new((/nfrq2/),float)
calc_option = 1

if(calc_option.eq.0) then 
    frq2 = fspan(0.001,0.5,nfrq2)
    var2 = linint1_n_Wrap(frq,var,False,frq2,0,0)
else
    do nfr = 0, nfrq2-1, 1
        n1 = nsta - 1 + nspc * nfr
        n2 = nsta - 1 + nspc *(nfr+1)
        n3 = (n1+n2)/2
        var2(nfr,:) = dim_avg_n_Wrap(var(n1:n2,:),0)
        frq2(nfr) = frq(n3)
    end do
    var2!0   = "frq"
    var2&frq = frq2
end if
printVarSummary(var2)
print(frq2(0:20))
print(frq(0:1000:100))
print(var2(0:20,5))
print(var(0:1000:100,5))

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
    res@vpHeightF        = 0.8
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = True
    res@tiYAxisOn       = True
    res@tiXAxisFont     = font_index 
    res@tiYAxisFont     = font_index 
    res@tiYAxisString   = "Frequency"
    res@tiXAxisString   = "Longitude(E)"

    res@tmYLMode        = "Explicit"
    res@tmYLValues      = fspan(0.0,0.5,6)
    res@tmYLLabels      = fspan(0.0,0.5,6)

    res@tmXBMode        = "Manual"
    res@tmXBTickStartF  = lonl(na)
    res@tmXBTickEndF    = lonr(na)
    res@tmXBTickSpacingF= lon_spcing 

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
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Vertical";"Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2
    res@pmLabelBarWidthF  = 0.1
    res@pmLabelBarHeightF = 0.8

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    if(area_option) then 
    res@gsnRightString  = lev(nl) + "hPa " + var_name(nv) +" ("+ lats(na)+"N)"   
    else
    res@gsnRightString  = lev(nl) + "hPa " + var_name(nv) +" ("+ lats(na)+" - "+latn(na)+"N)"   
    end if

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

var2 = var2/conform(var2,dim_avg_n(var2,1),0)-1
res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
    minlvl := (/-0.7/) ;used for East Asian moist transiation
    maxlvl := (/ 0.7/) ;used for CTRL-NUDG
    spacig := (/ 0.1/)
    
    ;minlvl := (/10 /) ;used for East Asian moist transiation
    ;maxlvl := (/490/) ;used for CTRL-NUDG
    ;spacig := (/30 /)

ip   = ind(frq.ge.0.001)    ; all indices for "short" periods

    ;res@cnFillPalette := "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    res@gsnLeftString   = pre_case(0) + ") " + case
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF      := spacig        ; set contour interval
    res@cnMaxLevelValF       := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF       := minlvl       ; set min contour level
    ;plot  = gsn_csm_hov(wks,var(ip,:),res)
    plot(nl)  = gsn_csm_hov(wks,var2,res)
    print(max(var2))
end do    
    ;resp@txString      = fig_title
    resp@txFont        = font_index
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    gsn_panel(wks,plot,(/2,2/),resp)
end

