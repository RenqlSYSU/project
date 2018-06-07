;*********************************************************************************
;filename:      draw_1p4x2_JJAvort_EOF.ncl
;input_data:    JJApreci_EOFn3.nc
;Draw the first three EOF spatial patterns and principal components of Asia JJA preci data
;every panel have three figures of spatial patterns (CTRL,NUDG,DIFF) 
;and one figure of time series which has two curves (CTRL,NUDG).
;                                            by Ql_Ren
;                                           2017/12/18
;*********************************************************************************
begin
ne   = 2
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-JJAvort_daily_EEOF_avg.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/1p8x3_JJAvort_daily_EEOF" + (ne+1) + "_avg"
case     = (/"CTRL","NG24h","NG6h"/)
pre_case = (/"(a","(b","(c"/)
fig_out  = "pdf" ;"eps";"x11"
unit     = "s~S~-1~N~"
var_name = (/"eof","eof_ts"/)

;plot area 
lats = 25 
latn = 40 
lonl = 70 
lonr = 130
lev  = 500

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
ndelay = getfilevardimsizes(f,"delay")
neval  = getfilevardimsizes(f,"neval")
ntime  = getfilevardimsizes(f,"time")
nlat   = getfilevardimsizes(f,"lat")
nlon   = getfilevardimsizes(f,"lon")
map_var  = new((/ncase,ndelay,nlat,nlon/),float)
ts_var   = new((/ncase,ndelay,ntime/),float)
pcvar    = new((/ncase/),float)
eval     = new((/ncase/),float)
ts_mean  = new((/ncase/),float)
sig      = new((/ncase/),logical)
plot     = new((/ncase*ndelay/),graphic)
line     = new((/ncase*ndelay/),graphic)
gtop     = new((/ncase*ndelay/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.04
   subtitle_font  = 0.04
   title_font     = 0.04
   label_font     = 0.04
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;================================================================
;read data
;====================================================================
do nf = 0,ncase-1,1
    f = addfile(filename(nf),"r")
    map_var(nf,:,:,:) = f->$var_name(0)$(ne,:,:,:)       ;(ncase,neval,nlat,nlon)
    pcvar(nf)   = map_var@pcvar(ne)
    eval(nf)    = map_var@eval(ne)
    sig(nf) = eofunc_north(pcvar(nf),ntime,True)
end do
subtitle_color = where(sig,"red","black") 
printVarSummary(map_var)

;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "landscape"

;------------set the attribute of the map----------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 10
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

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
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -0.04       ; set min contour level
    res@cnMaxLevelValF        = 0.04      ; set max contour level
    res@cnLevelSpacingF       = 0.005        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200

    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-2
    
    ;res@pmLabelBarZone               = 0
    ;res@pmLabelBarOrthogonalPosF     = 0.15
    ;res@pmLabelBarParallelPosF       = 0.5

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
    resp  = True    
    resp@gsnMaximize   = True
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    ;resp@gsnPanelRowSpec  = True
    resp@gsnPaperOrientation = "portrait";"landscape"
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Horizontal"
    resp@lbLabelStride    = 1
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@pmLabelBarWidthF  = 0.5
    resp@pmLabelBarHeightF = 0.05
    resp@lbLabelFontHeightF = tick_font-0.02
    resp@lbLabelFontThicknessF = tick_thick
    resp@lbBoxLineThicknessF  = mp_thick-1

do nd = 0,ndelay-1,1
    nn = nd + 1
    do nc = 0,ncase-1,1
        np = ncase*nd + nc
        res@gsnLeftString    = pre_case(nc) + nn + ") " + case(nc)
        res@gsnRightStringFontColor   = subtitle_color(nc)
        res@gsnRightString            = eval(nc)+"  "+pcvar(nc)+"%" 
        plot(np) = gsn_csm_contour_map_ce(wks,map_var(nc,nd,:,:),res)
        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(map_var(nc,nd,{lats:latn},{lonl:lonr})))
        print(min(map_var(nc,nd,{lats:latn},{lonl:lonr})))
    end do
end do    
gsn_panel(wks,plot,(/ndelay,ncase/),resp)
end

