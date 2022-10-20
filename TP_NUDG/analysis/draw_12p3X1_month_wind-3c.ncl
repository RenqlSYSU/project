;*********************************************************************
;filename:      draw_12p3X1_month_preci-3c_shad.ncl
;input_data:    CTRL-Clim_month_ave_preci.nc 
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2018/01/18
;Modification:
;20180308 modified for draw temperature or omega or q, but before draw you should 
;         modify some parameters: draw_var,  unit, var_name, months, nl, prob of read data
;         minlevel, max level and level spacing.
;*********************************************************************************
begin
draw_var = (/"U","V","wind"/)  ;this is the variable used to draw, can used as center string
unit     = "m/s" ;"K"   ;this can be used as right string 
;unit     = "m~S~2~N~/s~S~2~N~"   ;this can be used as right string 
;var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;months   = (/"DJF","MAM","JJA","SON"/)
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
lev  = (/850,500,200/)  ;can used as right string
nl = 1;1
nv = 2

;the path of the file
;path     = "/users/yangsong3/renql/project/TP_NUDG"
path    = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/analysis/fig/12p3X1_" + draw_var(nv) + "_" + var_name(0) + "-3c-" + lev(nl)
;fig_title  = "Seasonly " + draw_var(nv) + " Variance in " ;+season
fig_title  = "Climate Seasonal " + draw_var(nv) + " in " ;+season
    ;case   = (/"AMIP_CTRL","(CTRL-NUDG24h)/CTRL%","(CTRL-NUDG6h)/CTRL%"/)
    case   = (/"CTRL","NG6h-CTRL","NG24h-CTRL"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = -90 ;-15 ;-20
latn =  90 ;55  ;60
lonl =  0  ;30  ;20
lonr = 360 ;150 ;220
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
f = addfile(filename(0)+ draw_var(0) + ".nc","r")
ncase  = dimsizes(case)
nmonth = dimsizes(months)
nlat   = getfilevardimsizes(f,"lat")
nlon   = getfilevardimsizes(f,"lon")
nlev   = getfilevardimsizes(f,"lev")
var_u  = new((/ncase  ,nmonth,nlat,nlon/),float)
var_v  = new((/ncase  ,nmonth,nlat,nlon/),float)
vari_u = new((/ncase  ,nmonth,nlat,nlon/),float)
vari_v = new((/ncase  ,nmonth,nlat,nlon/),float)
prob_u = new((/ncase-1,nmonth,nlat,nlon/),float)
prob_v = new((/ncase-1,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/3/),graphic)
    line     = new((/3/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   vcref_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   vc_thick       = 1.5

;================================================================
;read data
;====================================================================
nv = 0
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf)+ draw_var(nv) + ".nc","r")
    var_u(nf,:,:,:)  = f->$var_name(0)$(:,{lev(nl)},:,:) ;get an array (ncase,nmonth,nlat,nlon)
    vari_u(nf,:,:,:) = f->$var_name(1)$(:,{lev(nl)},:,:);get an array (ncase,nmonth,nlat,nlon)
end do

nv = 1
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf)+ draw_var(nv) + ".nc","r")
    var_v(nf,:,:,:)  = f->$var_name(0)$(:,{lev(nl)},:,:) ;get an array (ncase,nmonth,nlat,nlon)
    vari_v(nf,:,:,:) = f->$var_name(1)$(:,{lev(nl)},:,:);get an array (ncase,nmonth,nlat,nlon)
end do

do nf = 0 ,ncase-2,1
;prob(nf,:,:,:) = ftest(var(0,:,:,:),27,var(nf+1,:,:,:),27,0)
prob_u(nf,:,:,:) = ttest(var_u(0,:,:,:),vari_u(0,:,:,:),27,var_u(nf+1,:,:,:),vari_u(nf+1,:,:,:),27,True,False)
prob_v(nf,:,:,:) = ttest(var_v(0,:,:,:),vari_v(0,:,:,:),27,var_v(nf+1,:,:,:),vari_v(nf+1,:,:,:),27,True,False)
var_u(nf+1,:,:,:)= var_u(nf+1,:,:,:) - var_u(0,:,:,:);/abs(var_u_temp)*100 ;get the difference between CTRL and NUDG
var_v(nf+1,:,:,:)= var_v(nf+1,:,:,:) - var_v(0,:,:,:);/abs(var_v_temp)*100 ;get the difference between CTRL and NUDG
end do

var_u(1:2,:,:,:) = mask(var_u(1:2,:,:,:),((prob_u.lt.siglvl).or.(prob_v.lt.siglvl)),True)
var_v(1:2,:,:,:) = mask(var_v(1:2,:,:,:),((prob_u.lt.siglvl).or.(prob_v.lt.siglvl)),True)
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
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@gsnScalarContour = False
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    res@cnMonoFillPattern= True;False
    res@cnFillPattern    = 0
    res@cnMonoFillColor  = True
    res@cnFillColor      = "gray"
    res@lbLabelBarOn     = False;True   ; turn off the label bar
   
;----------set the attribute of the vector-------------------
     res@vcGlyphStyle               = "CurlyVector"
     res@vcVectorDrawOrder          = "PostDraw"
     res@vcLineArrowThicknessF      = vc_thick
     res@vcMinDistanceF             = 0.02
     ;res@vcRefMagnitudeF            = 30.0             ; define vector ref mag
     res@vcRefLengthF               = 0.05           ; define length of vec ref
     res@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
     ;res@vcRefAnnoSide              = "Bottom"          
     ;res@vcRefAnnoJust              = "TopRight"
     res@vcRefAnnoOrthogonalPosF    = 0.5           ; move ref vector along Y-axis
     res@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
     res@vcRefAnnoPerimOn           = False         ;close the border
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
    res_gtop@cnLineThicknessF = mp_thick
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
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font

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

    vcref         = (/(/10,3,3/), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/10,3,3/), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/25,3,3/)/)   ;the level of 200hpa    ctrl,nudg,diff
do nm = 0 ,nmonth-1,1
    do nc = 0, ncase-1,1
        res@gsnLeftString   = case(nc)
        ;if(nc.gt.0) then
        ;    res@gsnScalarContour = True
        ;    res@vcRefMagnitudeF       = vcref(nl,nc)                ; define vector ref mag
        ;    res@vcRefAnnoString2      = vcref(nl,nc)+"m/s"
        ;    res@gsnRightString  = "siglvl: " + siglvl
        ;    plot(nc)     = gsn_csm_vector_scalar_map_ce(wks,var_u(nc,nm,:,:),var_v(nc,nm,:,:),prob(nc-1,nm,:,:),res)
        ;else
            res@gsnScalarContour = False
            res@vcRefMagnitudeF       = vcref(nl,nc)             ; define vector ref mag
            res@vcRefAnnoString2      = vcref(nl,nc)+"m/s"
            res@gsnRightString   = unit + " " + lev(nl) + "hPa" ;draw_var(nv) 
            plot(nc)     = gsn_csm_vector_map_ce(wks,var_u(nc,nm,:,:),var_v(nc,nm,:,:),res)
        ;end if
        line(nc) = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
        print(max(var_u(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(var_u(nc,nm,{lats:latn},{lonl:lonr})))
    end do 

    resp@txString      = fig_title + months(nm) 
    gsn_panel(wks,plot,(/3,1/),resp)
end do ;for the month 
end

