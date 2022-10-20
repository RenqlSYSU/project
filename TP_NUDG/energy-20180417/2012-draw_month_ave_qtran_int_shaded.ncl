;*********************************************************************
;filename:      draw_1p2x3_qtran_int.ncl
;input_data:    CTRL-JJA_ano_850hpa_qu.nc 
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
draw_var  = (/"Total","Steady","Transient"/) 
var_name  = (/"qu_ave","qv_ave","div_ave"/)
vari_name = (/"qu_vari","qv_vari","div_vari"/)
unit      = "kg/(s*m)" ;"K"   ;this can be used as right string 
diffterm  = False;True ;
sig_option= False;True ;if False, only the significant values are shaded 
barshare  = True ;False;
barvert   = True ;False;

months   = (/"MJ","JA"/)
imonth = ispan(0,1,1);(/4,5,6,7/)
;months   = (/"JF","MA","MJ","JA"/)
;imonth = ispan(2,3,1);(/4,5,6,7/)
;months  = (/"May","Jun","Jul","Aug"/)
nmonth = dimsizes(imonth)
month_s = 5
month_e = 8

;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG6h-CTRL"/) ;,"NG24h-CTRL"
    gtopfile= path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
end if
if(case_option.eq.1) then
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"CTRL-noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
end if
if(case_option.eq.2) then
    case1   = (/"F2000","NGTR"/)       ;NGTR    
    case2   = (/"F2000","NGTR-F2000"/) ;NGTR
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 30
    if(diffterm) then
    case1   := (/"Ano2"/)       ;NGTR    
    case2   := (/"NGTR-F2000"/) ;NGTR
    draw_var   := (/"q'V","qV'"/) 
    end if
end if
filename = path + "/energy-20180417/mdata/" + case1 + "-Clim_2month_ave_qtran_int.nc" ;,"NUDG24h"
fig_name = path + "/energy-20180417/fig/" + case1(0) + "-qtran_int_2month"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p"/)
;fig_out    ="pdf";"x11"
fig_out    ="eps";"x11"

;area: east Asia
lats = (/10 ,-15/) ;-20
latn = (/55 ,55 /) ;60
lonl = (/90 ,30 /) ;20
lonr = (/150,150/) ;220
siglvl = 0.05
na = 0

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case1)
nvar1  = dimsizes(draw_var)
nvar2  = dimsizes(var_name)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)

;the number of the panel plot
    plot     = new((/5*nmonth/),graphic)
    plotvc   = new((/5*nmonth/),graphic)
    wind2    = new((/5*nmonth/),graphic)
    line     = new((/5*nmonth/),graphic)
    gtop     = new((/5*nmonth/),graphic)

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
var  = new((/ncase,nvar2,nvar1,nmonth,nlat,nlon/),double)
vari = new((/ncase,nvar2,nvar1,nmonth,nlat,nlon/),double)

do nf = 0 ,ncase-1, 1
f = addfile(filename(nf),"r")
do nv = 0,nvar2-1,1
    print("read " + case1(nf) + " " + var_name(nv))
    var(nf,nv,:,:,:,:)   =  f->$var_name(nv)$(:,imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    if(.not.diffterm) then 
    vari(nf,nv,:,:,:,:)  = f->$vari_name(nv)$(:,imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end do
end do

if(diffterm) then
var(:,2,:,:,:,:) = var(:,2,:,:,:,:)*100000
else
prob = new((/ncase-1,nvar2,nvar1,nmonth,nlat,nlon/),double)
do nf = 0 ,ncase-2,1
prob(nf,:,:,:,:,:) = ttest(var(0,:,:,:,:,:),vari(0,:,:,:,:,:),nsample,var(nf+1,:,:,:,:,:),vari(nf+1,:,:,:,:,:),nsample,True,False)
end do

var(:,2,:,:,:,:) = var(:,2,:,:,:,:)*100000
;var(:,2,:,:,:,:) = -uv2dv_cfd ( var(:,0,:,:,:,:), var(:,1,:,:,:,:), var&lat, var&lon, 2)*100000
do nc = 1, ncase-1, 1
    var(nc,:,:,:,:,:) =  var(nc,:,:,:,:,:) - var(0,:,:,:,:,:) 
end do
printVarSummary(var)

uvar = mask(var(1:(ncase-1),0,:,:,:,:),((prob(:,0,:,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:,:).lt.siglvl)),False) ;not significance
vvar = mask(var(1:(ncase-1),1,:,:,:,:),((prob(:,0,:,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:,:).lt.siglvl)),False)
var(1:(ncase-1),0,:,:,:,:) = mask(var(1:(ncase-1),0,:,:,:,:),((prob(:,0,:,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:,:).lt.siglvl)),True)
var(1:(ncase-1),1,:,:,:,:) = mask(var(1:(ncase-1),1,:,:,:,:),((prob(:,0,:,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:,:).lt.siglvl)),True)
copy_VarMeta(var(0,0,0,0,:,:),uvar(0,0,0,:,:))
copy_VarMeta(var(0,0,0,0,:,:),vvar(0,0,0,:,:))
if(sig_option) then 
prob(:,2,:,:,:,:) = mask(prob(:,2,:,:,:,:),prob(:,2,:,:,:,:).lt.siglvl,True)
else
var(1:(ncase-1),2,:,:,:,:) = mask(var(1:(ncase-1),2,:,:,:,:),prob(:,2,:,:,:,:).lt.siglvl,True)
;var(1:(ncase-1),2,:,:,:,:) = mask(var(1:(ncase-1),2,:,:,:,:),((prob(:,0,:,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:,:).lt.siglvl).or.(prob(:,2,:,:,:,:).lt.siglvl)),True)
end if
end if

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   vc_thick       = 3.0
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
;    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = 15
    res@gsnMinorLatSpacing = 15
    res@gsnMajorLonSpacing = 20
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
    ;res@tmXBLabelsOn         = False              ; do not draw bottom labels
    ;res@tmXBOn               = False              ; no bottom tickmarks
    ;res@tmYLLabelsOn         = False              ; no right labels
    ;res@tmYLOn               = False              ; no right tickmarks
    ;res@tmXTLabelsOn         = False              ; do not draw bottom labels
    ;res@tmXTOn               = False              ; no bottom tickmarks
    ;res@tmYLLabelsOn         = False              ; no right labels
    ;res@tmYLOn               = False              ; no right tickmarks
    
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
    res@gsnLeftString      = ""
    res@gsnRightString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    
;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
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
    
;----------set the attribute of the vector-------------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMinDistanceF             = 0.028
    resvc@vcLineArrowHeadMinSizeF    = 0.013
    resvc@vcLineArrowHeadMaxSizeF    = 0.013
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoOn                = False
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.3           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5           ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    ;resvc@vcRefAnnoString2           = "30m/s"

    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor           = "grey"

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
    res_gtop@cnLevels         = (/3000/)
    res_gtop@gsnCenterString      = ""
    res_gtop@gsnLeftString    = " "
    res_gtop@gsnRightString   = " "

;----------Set res sig-----------
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

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"  ;"portrait"
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
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
    res@lbBoxLineThicknessF  = mp_thick-2
    res@pmLabelBarZone       = 0
    if(barvert) then
    res@lbOrientation        = "Vertical";
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.17
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.08
    res@pmLabelBarOrthogonalPosF     = 0.15 ;y direction
    res@pmLabelBarParallelPosF       = 0.5
    end if
    end if
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(diffterm) then 
    vcref  = (/(/50,50,20/),(/50,50,20/),(/50,50,20/)/)  ;the level of 850hpa    ctrl,nudg,diff
else
    vcref  = (/(/200,200,50/),(/50,50,20/),(/50,50,20/)/)  ;the level of 850hpa    ctrl,nudg,diff
end if

    ;minlvl = (/(/15 ,45 , 0/),(/-80,-16,-16/),(/-80,-16,-16/)/) ;used for East Asian Shaded ration
    ;maxlvl = (/(/615,120,30/),(/ 80, 16, 16/),(/ 80, 16, 16/)/)
    ;spacig = (/(/40 ,5  , 2/),(/ 10, 2 , 2 /),(/ 10, 2 , 2 /)/)
    res@cnFillPalette = "BlueDarkRed18"
    res@cnFillColors  = (/0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17/) ; If the color indexes are derived from cnFillPalette the first color comes from element 0
    maxlvl = (/(/ 7, 7, 7/),(/ 7, 7, 7/),(/ 7, 7, 7/)/) ;used for East Asian moist transiation
    minlvl = (/(/-7,-7,-7/),(/-7,-7,-7/),(/-7,-7,-7/)/) ;used for CTRL-NUDG
    spacig = (/(/ 1, 1, 1/),(/ 1, 1, 1/),(/ 1, 1, 1/)/)

np = -1
do nc = 0, ncase-1,1
do nv = 1 ,nvar1-1,1
do nm = 0, nmonth-1,1
        np = np + 1
        ;np = nvar*nc + nv
        ;np = nmonth*nc + nm
        ;np = nmonth*nv + nm
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(imonth(nm));+ " " + draw_var(nv)
        res@gsnRightString  = draw_var(nv)
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF      := spacig(nc,nv)        ; set contour interval
        res@cnMaxLevelValF       := maxlvl(nc,nv) ;32    ; set max contour level
        res@cnMinLevelValF       := minlvl(nc,nv)       ; set min contour level
        resvc@vcRefAnnoOn       = True
        resvc@vcRefMagnitudeF   = vcref(nc,nv)                ; define vector ref mag
        resvc@vcRefAnnoString2  = vcref(nc,nv); + " " + unit
        resvc2@vcRefMagnitudeF  = vcref(nc,nv)                ; define vector ref mag
        resvc2@vcRefAnnoString2 = vcref(nc,nv) + " " + unit
        plot(np)   = gsn_csm_contour_map_ce(wks,var(nc,2,nv,nm,:,:),res) 
        if(.not.diffterm) then 
        if(sig_option.and.nc.gt.0) then 
        sig(np) = gsn_csm_contour(wks,prob(nc-1,2,nv,nm,:,:),res_sig)
        overlay(plot(np),sig(np))
        end if
        
        if(nc.gt.0)
        wind2(np)  = gsn_csm_vector(wks,uvar(nc-1,nv,nm,:,:),vvar(nc-1,nv,nm,:,:),resvc2)
        overlay(plot(np),wind2(np))
        end if
        end if
        plotvc(np) = gsn_csm_vector(wks,var(nc,0,nv,nm,:,:),var(nc,1,nv,nm,:,:),resvc)
        overlay(plot(np),plotvc(np))
        
        gtop(np)   = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        ;line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var(nc,2,nv,nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        print(min(var(nc,2,nv,nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
    end do 
end do 
    ;res1 = True
    ;res2 = True
    ;res1@gsnMaximize = True
    ;attach1 = gsn_attach_plots(plot(0),plot(1:3),res1,res2)
    ;attach2 = gsn_attach_plots(plot(4),plot(5:7),res1,res2)
    ;gsn_panel(wks,plot(0:7:4),(/ncase,1/),resp)
    ;gsn_panel(wks,plot,(/ncase,nvar/),resp)
 ;   gsn_panel(wks,plot,(/ncase,nmonth/),resp)
end do
;    gsn_panel(wks,plot,(/nvar1,nmonth/),resp)
    gsn_panel(wks,plot,(/4,2/),resp)
end 

