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

if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG6h"/) ;,"NG24h"
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

area = (/"EA","TP","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 25 ,25  , 22  , 30  /)
latn = (/ 40 ,40  , 30  , 45  /)
lonl = (/ 30 ,70  ,107  ,115  /)
lonr = (/150 ,105 ,130  ,130  /)
lev  = (/850,500,200/)
var_name = (/"U","V","T"/)
nl = 2
na = 0
nv = 1
filein   = path + "/mdata/" + case1 + "-allyear_spec_" + var_name(nv) + "-" + area(na) + ".nc"
fig_name = path + "/fig/" + case1(0) + "-allyear_spec_" + var_name(nv) + "-" + area(na) + lev(nl)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

f    = addfile(filein(0),"r")
spcx = f->spcx(:,{lev(nl)},:)
frq  = f->frq
nfrq = dimsizes(spcx&frq)
nlon = dimsizes(spcx&lon)
ncase = dimsizes(case1)
var  = new((/ncase,nfrq,nlon/),float)
plot = new((/ncase/),graphic)

do nc = 0, ncase-1, 1
    f    = addfile(filein(nc),"r")
    var(nc,:,:) = f->spcx(:,{lev(nl)},:)
end do
printVarSummary(var)
copy_VarMeta(spcx,var(0,:,:))
    
;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22
   
   wks = gsn_open_wks(fig_out,fig_name)  

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    res@vpHeightF        = 0.25
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiYAxisOn       = False
    res@tiYAxisString   = "Frequency"
    res@tiXAxisString   = "Longitude"

;    res@tmYLMode        = "Explicit"
;    res@tmYLValues      = fspan(0.0,0.5,0.1)
;    res@tmYLLabels      = fspan(0.0,0.5,0.1)

    ;res@tmYLMode        = "Manual"
    ;res@tmYLTickStartF  = lats(na)
    ;res@tmYLTickEndF    = latn(na)
    ;res@tmYLTickSpacingF= lat_spacing

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
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 0         ; set min contour level
    res@cnMaxLevelValF        = 160       ; set max contour level
    res@cnLevelSpacingF       = 10        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString  = var_name(nv) +" "+ lats(na)+" - "+latn(na)+"N"   

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

;    minlvl = (/0 ,0 ,-4 /) ;used for East Asian moist transiation
;    maxlvl = (/16,16, 4 /) ;used for CTRL-NUDG
;    spacig = (/1 ,1 ,0.5/)
    
ip   = ind(frq.ge.0.001)    ; all indices for "short" periods

do nc = 0, ncase-1,1
if(nc.le.2) then
    res@cnFillPalette := "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
else
    res@cnFillPalette := "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
end if 
    res@gsnLeftString   = pre_case(nc) + ") " + case1(nc)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnLevelSpacingF      := spacig(nc)        ; set contour interval
    ;res@cnMaxLevelValF       := maxlvl(nc) ;32    ; set max contour level
    ;res@cnMinLevelValF       := minlvl(nc)       ; set min contour level
    plot(nc)  = gsn_csm_hov(wks,var(nc,ip,:),res)
    print(max(var(nc,ip,:)))
    print(min(var(nc,ip,:)))
end do 

    gsn_panel(wks,plot,(/ncase,1/),resp)
end

