;*********************************************************************
;filename:      draw_1p2x4_qtran_ano_850hpa.ncl
;input_data:    CTRL-JJA_ano_850hpa_qu.nc 
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
file_pos = (/"qu","qv"/)  ;this is the variable used to draw, can used as center string
draw_var = (/"(qu)~s~-~n~&(qv)~s~-~n~","q~s~-~n~u~s~-~n~&q~s~-~n~v~s~-~n~","10~s~-9~n~ q~s~-~n~u'&q~s~-~n~v'","~s~-9~n~  q'u~s~-~n~&q'v~s~-~n~"/) 
var_name = (/"ave_qu","aveq_aveu","aveq_u","q_aveu"/);,(/"ave_qv","aveq_avev","aveq_v","q_avev"/)/)
unit     = "g/(s*Pa*m)" ;"K"   ;this can be used as right string 

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h"/) + "-JJA_ano_850hpa_"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/1p2x4_qtran_JJAano_850hpa_"
    case   = (/" CTRL"," NG6h-CTRL"/)
pre_case   = (/"(a","(b"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
lev  = (/850/)
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case)
nvar   = dimsizes(draw_var)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var_u  = new((/ncase,nvar,nlat,nlon/),float)
var_v  = new((/ncase,nvar,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase*nvar/),graphic)
    line     = new((/ncase*nvar/),graphic)
    gtop     = new((/ncase*nvar/),graphic)

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
   vcref_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   vc_thick       = 1.5
   font_index     = 22

;================================================================
;read data
;====================================================================
npos = 0
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf)+ file_pos(npos) + ".nc","r")
    do nv = 0,nvar-1,1
    var_u(nf,nv,:,:)  = f->$var_name(nv)$(:,:) ;get an array (ncase,nmonth,nlat,nlon)
    end do
end do

npos = 1
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf)+ file_pos(npos) + ".nc","r")
    do nv = 0,nvar-1,1
    var_v(nf,nv,:,:)  = f->$var_name(nv)$(:,:) ;get an array (ncase,nmonth,nlat,nlon)
    end do
end do

var_u(1,:,:,:) = var_u(1,:,:,:) - var_u(0,:,:,:)
var_v(1,:,:,:) = var_v(1,:,:,:) - var_v(0,:,:,:)
var_u(:,2:3,:,:) = var_u(:,2:3,:,:)*1000000000
var_v(:,2:3,:,:) = var_v(:,2:3,:,:)*1000000000
printVarSummary(var_u)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
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

;----------set the attribute of the vector-------------------
     res@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
     res@vcVectorDrawOrder          = "PostDraw"
     res@vcLineArrowThicknessF      = vc_thick
     res@vcMinDistanceF             = 0.03
     ;res@vcRefMagnitudeF            = 30.0             ; define vector ref mag
     res@vcRefLengthF               = 0.05           ; define length of vec ref
     res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
     ;res@vcRefAnnoSide              = "Bottom"          
     ;res@vcRefAnnoJust              = "TopRight"
     res@vcRefAnnoOrthogonalPosF    = -0.67           ; move ref vector along Y-axis
     res@vcRefAnnoParallelPosF      = 0            ; move ref vector along X-axis
     res@vcRefAnnoPerimOn           = False         ;close the border
     res@vcRefAnnoFont              = font_index
     res@vcRefAnnoFontHeightF       = vcref_font
     res@vcRefAnnoString1On         = False
     res@vcRefAnnoString2On         = True
     ;res@vcRefAnnoString2           = "30m/s"

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
    res_gtop@cnLevels         = (/15000/)
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
    res@gsnStringFontHeightF  = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait";"landscape"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    vcref         = (/(/10,10,10,1/),(/1,1,10,1/)/)  ;the level of 850hpa    ctrl,nudg,diff

do nv = 0 ,nvar-1,1
    do nc = 0, ncase-1,1
        nn = nc +1
        np = ncase*nv + nc
        res@gsnLeftString   = pre_case(nc) + nn + case(nc)
        res@gsnRightString  = lev + "hPa " + draw_var(nv)
        res@vcRefMagnitudeF       = vcref(nc,nv)                ; define vector ref mag
        res@vcRefAnnoString2      = vcref(nc,nv) + " " + unit
        plot(np) = gsn_csm_vector_map_ce(wks,var_u(nc,nv,:,:),var_v(nc,nv,:,:),res)
        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var_u(nc,nv,{lats:latn},{lonl:lonr})))
        print(min(var_u(nc,nv,{lats:latn},{lonl:lonr})))
    end do 
end do 
    gsn_panel(wks,plot,(/nvar,ncase/),resp)
end

