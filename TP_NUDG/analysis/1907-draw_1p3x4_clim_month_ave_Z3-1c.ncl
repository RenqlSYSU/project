;*********************************************************************
;Draw monthly preci in every month(3*4) of one case, just simply view 
;the distribution of a variable.
;                                                     by Ql Ren
;                                                    2019/07/12
;*********************************************************************************
begin
draw_var = (/"preci" ,"EKE"        ,"U"  ,"Z3","PS"           ,"U10"/)  ;this is the variable used to draw, can used as center string
unit     = (/"mm/day","10~S~4~N~ J","m/s","m" ,"hPa(+1000hPa)","m/s"/)    ;this can be used as right string 
;var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;months   = (/"DJF","MAM","JJA","SON"/)
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
lev = (/925,850,500,200/)
nl = 3
nv = 2
case_option = 0

;the path of the file
path    = "/home/ys17-19/renql/project/TP_NUDG"
if(case_option.eq.1) then 
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
else
    case1   = (/"CTRL","NUDG6h","NUDG24h"/)
    case2   = (/" CTRL"," NG6h-CTRL", " NG24h-CTRL"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
end if

nc = 0
filename  = path + "/analysis/mdata/" + case1(nc) + "-Clim_month_ave_" + draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/1p3x4_" + case2(nc) + "_" + lev(nl) + draw_var(nv) + "_" + var_name(0) + "-2c" 
;fig_title  = "Seasonal " + draw_var(nv) + " variance in " ;+season
fig_title  = "Climate " + draw_var(nv) + " (" + unit(nv) + ")"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = -15 ;20 ;-20
latn = 55  ;55 ;60
lonl = 30  ;30 ;20
lonr = 150 ;150;220

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nmonth = dimsizes(months)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var  = new((/nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/nmonth/),graphic)
    gtop     = new((/nmonth/),graphic)
    line     = new((/nmonth/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
if(nv.ge.4) then 
    var  = f->$var_name(0)$(:,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
else
    var  = f->$var_name(0)$(:,{lev(nl)},{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
end if

if(nv.eq.1) then
    var = var/10000
end if

if(nv.eq.4) then
    var = var/100-1000
end if

printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
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
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    ;res@cnLineThicknessF = 0
    ;res@gsnContourZeroLineThicknessF = 4.0
    ;res@cnLineColor      = "blue"
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 1          ; set min contour level
    ;res@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5
;----------Set res sig-----------
    res_sig                  = True
    res_sig@gsnDraw          = False
    res_sig@gsnFrame         = False
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnInfoLabelOn    = False
    ;res_sig@cnLineThicknessF = 2
    ;res_sig@cnLineColor     = "black"
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

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 4.0
    ;lnres@gsLineLabelString    = "TP"

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = 4.0
    res_gtop@cnLineColor      = "black"
    res_gtop@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_gtop@cnLevels         = (/1500/)
    res_gtop@gsnCenterString      = ""
    res_gtop@gsnLeftString    = " "
    res_gtop@gsnRightString   = " "

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(nv.eq.0) then ;preci
    minlvl = (/ 1 /);-8,-8,-8,-8
    maxlvl = (/16 /);8 ,8 ,8 ,8 
    spacig = (/ 1 /);1 ,1 ,1 ,1 
end if

if(nv.eq.1) then ;EKE
    minlvl = (/ 2 /);-8,-8,-8,-8
    maxlvl = (/34 /);8 ,8 ,8 ,8 
    spacig = (/ 2 /);1 ,1 ,1 ,1 
end if

if(nv.eq.2) then ;U
    minlvl = (/-4 /); 2 -8,-8,-8,-8
    maxlvl = (/ 4 /);50 8 ,8 ,8 ,8 
    spacig = (/0.5/); 3 1 ,1 ,1 ,1 
end if

if(nv.eq.3) then ;Z3
    minlvl = (/5600 /);-8,-8,-8,-8
    maxlvl = (/5920 /);8 ,8 ,8 ,8 
    spacig = (/ 20  /);1 ,1 ,1 ,1 
end if

if(nv.eq.4) then ;PS
    minlvl = (/-8/);-8,-8,-8,-8
    maxlvl = (/ 8/);8 ,8 ,8 ,8 
    spacig = (/ 1/);1 ,1 ,1 ,1 
end if

if(nv.eq.5) then ;U10
    minlvl = (/-4 /); 2 -8,-8,-8,-8
    maxlvl = (/ 4 /);50 8 ,8 ,8 ,8 
    spacig = (/0.5/); 3 1 ,1 ,1 ,1 
end if


if(nv.lt.1) then
   res@cnFillPalette = "precip3_16lev"
   ;blue and yellow and red, has 17 colors,2-18,2 is white
else
    res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
end if 
do nm = 0,nmonth-1,1
    res@gsnLeftString   = months(nm) 
    res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"
    res@gsnSpreadColors       = False              ; turn of the full color map
    if(nv.gt.6) then
    res@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
    else
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig        ; set contour interval
    res@cnMaxLevelValF        := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl       ; set min contour level
    end if
    plot(nm)     = gsn_csm_contour_map_ce(wks,var(nm,:,:),res)
    gtop(nm) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(nm),gtop(nm))
    ;line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
    print(max(var(nm,{lats:latn},{lonl:lonr})))
    print(min(var(nm,{lats:latn},{lonl:lonr})))
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/3,4/),resp)
end

