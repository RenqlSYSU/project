;*********************************************************************
;filename:      draw_1p2x3_qtran_int.ncl
;input_data:    CTRL-JJA_ano_850hpa_qu.nc 
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   vc_thick       = 2.5
   font_index     = 22

draw_var = (/"climate","steady","transient"/) 
var_name_u = (/"qu","ste_qu","edd_qu"/)
var_name_v = (/"qv","ste_qv","edd_qv"/)
vari_name  = (/"qu","qv","q"    ,"u"    ,"v"    ,"edd_qu","edd_qv"/) + "_vari"
ave_name   = (/"qu","qv","q_ave","u_ave","v_ave","edd_qu","edd_qv"/)
unit     = "kg/(s*m)" ;"K"   ;this can be used as right string 

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h"/) + "-Clim_JJAave_qtran_int.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/1p2x3_qtran_int_JJA"
    case   = (/" CTRL"," NG6h-CTRL"/)
pre_case   = (/"(a","(b"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
siglvl = 0.05
ntime  = 27*92

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case)
nvar   = dimsizes(draw_var)
nvari  = dimsizes(vari_name)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var_u  = new((/ncase,nvar ,nlat,nlon/),float)
var_v  = new((/ncase,nvar ,nlat,nlon/),float)
ave    = new((/ncase,nvari,nlat,nlon/),float) ;store the average of q,u,v
vari   = new((/ncase,nvari,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase*nvar/),graphic)
    plotvc   = new((/ncase*nvar/),graphic)
    line     = new((/ncase*nvar/),graphic)
    gtop     = new((/ncase*nvar/),graphic)

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
do nf = 0 ,ncase-1, 1
f = addfile(filename(nf),"r")
do nv = 0,nvar-1,1
    print("read " + case(nf) + " " + var_name_u(nv))
    print("read " + case(nf) + " " + var_name_v(nv))
    var_u(nf,nv,:,:)  = f->$var_name_u(nv)$({lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    var_v(nf,nv,:,:)  = f->$var_name_v(nv)$({lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
end do
do nv = 0,nvari-1,1
    print("read " + case(nf) + " " + vari_name(nv))
    print("read " + case(nf) + " " + ave_name(nv))
    vari(nf,nv,:,:)  = f->$vari_name(nv)$({lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
     ave(nf,nv,:,:)  = f->$ave_name(nv)$({lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
end do
end do

prob = ttest(ave(0,:,:,:),vari(0,:,:,:),ntime,ave(1,:,:,:),vari(1,:,:,:),ntime,True,False)
var_u(1,:,:,:) = var_u(1,:,:,:) - var_u(0,:,:,:)
var_v(1,:,:,:) = var_v(1,:,:,:) - var_v(0,:,:,:)
var_v(1,(/0,2/),:,:) = mask(var_v(1,(/0,2/),:,:),((prob((/0,5/),:,:).lt.siglvl).or.(prob((/1,6/),:,:).lt.siglvl)),True)
var_u(1,(/0,2/),:,:) = mask(var_u(1,(/0,2/),:,:),((prob((/0,5/),:,:).lt.siglvl).or.(prob((/1,6/),:,:).lt.siglvl)),True)
var_u(1,1,:,:) = mask(var_u(1,1,:,:),((prob(2,:,:).lt.siglvl).or.(prob(3,:,:).lt.siglvl).or.(prob(4,:,:).lt.siglvl)),True)
var_v(1,1,:,:) = mask(var_v(1,1,:,:),((prob(2,:,:).lt.siglvl).or.(prob(3,:,:).lt.siglvl).or.(prob(4,:,:).lt.siglvl)),True)

printVarSummary(var_u)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
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
    res@gsnMajorLatSpacing = 15
    res@gsnMinorLatSpacing = 15
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

    res@gsnCenterString = ""
    res@gsnLeftString   = ""
    res@gsnRightString  = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

;----------set the attribute of the vector-------------------
     res@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
     res@vcVectorDrawOrder          = "PostDraw"
     res@vcLineArrowThicknessF      = vc_thick
     res@vcMinDistanceF             = 0.035
     ;res@vcRefMagnitudeF            = 30.0             ; define vector ref mag
     res@vcRefLengthF               = 0.05           ; define length of vec ref
     res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
     ;res@vcRefAnnoSide              = "Bottom"          
     ;res@vcRefAnnoJust              = "TopRight"
     res@vcRefAnnoOrthogonalPosF    = -0.75           ; move ref vector along Y-axis
     res@vcRefAnnoParallelPosF      = 0.2            ; move ref vector along X-axis
     res@vcRefAnnoPerimOn           = False         ;close the border
     res@vcRefAnnoFont              = font_index
     res@vcRefAnnoFontHeightF       = vcref_font
     res@vcRefAnnoString1On         = False
     res@vcRefAnnoString2On         = True
     ;res@vcRefAnnoString2           = "30m/s"
     res@vcLineArrowHeadMinSizeF   = 0.015
     res@vcLineArrowHeadMaxSizeF   = 0.015

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
    resp@gsnPaperOrientation = "landscape"  ;"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    vcref  = (/(/200,200,20/),(/40,40,5/)/)  ;the level of 850hpa    ctrl,nudg,diff

do nc = 0, ncase-1,1
    do nv = 0 ,nvar-1,1
        nn = nv +1
        np = nvar*nc + nv
        res@gsnLeftString   = pre_case(nc) + nn + ")" + case(nc)
        res@gsnRightString  = draw_var(nv)
        res@vcRefMagnitudeF   = vcref(nc,nv)                ; define vector ref mag
        res@vcRefAnnoString2  = vcref(nc,nv) + " " + unit
        plot(np) = gsn_csm_vector_map_ce(wks,var_u(nc,nv,:,:),var_v(nc,nv,:,:),res)
        gtop(np)   = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var_u(nc,nv,{lats:latn},{lonl:lonr})))
        print(min(var_u(nc,nv,{lats:latn},{lonl:lonr})))
    end do 
end do 
    gsn_panel(wks,plot,(/ncase,nvar/),resp)
end

