;*********************************************************************
;filename:      draw_1p2x2_quv_TMQ_qdivg_ave_map-2c.ncl
;input_data:    CTRL-Clim_month_ave_*.nc 
;
;1 read JJA data from four files season_ave(2,:,:) and season_vari(2,:,:)
;2 do the significance test
;3 define the levels then draw the plot 
;                                                     by Ql Ren
;                                                    2018/03/28
;*********************************************************************************
begin
draw_var = (/"qu&qv&TMQ(kg/m~S~2~N~)","qdivg(10~S~-3~N~ g/(s*m*m))"/)  ;this is the variable used to draw, can used as center string
draw_varname = (/"TMQ","qdivg","qu","qv"/)
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months   = (/"DJF","MAM","JJA","SON"/)
nm = 2

;the path of the file
;path     = "/home/ys17-19/renql/project/TP_NUDG"
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/analysis/mdata/" + (/"CTRL","NUDG6h"/) + "-Clim_month_ave_" ;+ draw_varname(nv) + ".nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/phy_panel-20180328/fig/1p2x2_quv_TMQ_qdivg_" + months(nm) + var_name(0) 
    case   = (/" CTRL"," NG6h-CTRL"/)
pre_case   = (/"(a","(b"/)
fig_out    = "eps";"pdf";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
ncase  = dimsizes(case)
nmonth = dimsizes(months)
nvar   = dimsizes(draw_var)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var    = new((/ncase  ,nvar,nlat,nlon/),float)
var_u  = new((/ncase       ,nlat,nlon/),float)
var_v  = new((/ncase       ,nlat,nlon/),float)
vari   = new((/ncase  ,nvar,nlat,nlon/),float)
vari_u = new((/ncase       ,nlat,nlon/),float)
vari_v = new((/ncase       ,nlat,nlon/),float)
prob   = new((/ncase-1,nvar,nlat,nlon/),float)
prob_u = new((/ncase-1     ,nlat,nlon/),float)
prob_v = new((/ncase-1     ,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase*nvar/),graphic)
    line     = new((/ncase*nvar/),graphic)
    gtop     = new((/ncase*nvar/),graphic)
    plot2    = new((/(ncase-1)*nvar/),graphic)
    plot_uv   = new((/ncase/),graphic)
    plot_uv2  = new((/ncase-1/),graphic)

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
   vcref_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 3.0
   mp_thick       = 5.0
   font_index     = 22

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-1, 1
do nv = 0 ,nvar-1,1
    f = addfile(filename(nf)+draw_varname(nv)+".nc","r")
     var(nf,nv,:,:)  = f->$var_name(0)$(nm,:,:) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,nv,:,:)  = f->$var_name(1)$(nm,:,:) ;get an array (ncase,nmonth,nlat,nlon)
end do
f = addfile(filename(nf)+draw_varname(2)+".nc","r")
 var_u(nf,:,:)  = f->$var_name(0)$(nm,:,:) ;get an array (ncase,nmonth,nlat,nlon)
vari_u(nf,:,:)  = f->$var_name(1)$(nm,:,:) ;get an array (ncase,nmonth,nlat,nlon)

f = addfile(filename(nf)+draw_varname(3)+".nc","r")
 var_v(nf,:,:)  = f->$var_name(0)$(nm,:,:) ;get an array (ncase,nmonth,nlat,nlon)
vari_v(nf,:,:)  = f->$var_name(1)$(nm,:,:) ;get an array (ncase,nmonth,nlat,nlon)
end do

do nf = 0 ,ncase-2,1
prob(nf,:,:,:) = ttest(var(0,:,:,:),vari(0,:,:,:),27,var(nf+1,:,:,:),vari(nf+1,:,:,:),27,True,False)
prob_u(nf,:,:) = ttest(var_u(0,:,:),vari_u(0,:,:),27,var_u(nf+1,:,:),vari_u(nf+1,:,:),27,True,False)
prob_v(nf,:,:) = ttest(var_v(0,:,:),vari_v(0,:,:),27,var_v(nf+1,:,:),vari_v(nf+1,:,:),27,True,False)
var(nf+1,:,:,:)= var(nf+1,:,:,:) - var(0,:,:,:) ;get the difference between CTRL and NUDG
var_u(nf+1,:,:)= var_u(nf+1,:,:) - var_u(0,:,:) ;get the difference between CTRL and NUDG
var_v(nf+1,:,:)= var_v(nf+1,:,:) - var_v(0,:,:) ;get the difference between CTRL and NUDG
end do

var_u   = var_u/1000
var_v   = var_v/1000
var(:,1,:,:) = var(:,1,:,:)*1000     
;var(:,2,:,:) = var(:,2,:,:)/100      ;convert unit from m to 100 m

;var(1:2,:,:,:) = mask(var(1:2,:,:,:),prob.lt.siglvl,True)
prob   = mask(prob,prob.lt.siglvl,True)
u_test = var_u(1,:,:)
v_test = var_v(1,:,:)
u_test = mask(var_u(1,:,:),((prob_u(0,:,:).lt.siglvl).or.(prob_v(0,:,:).lt.siglvl)),True)
v_test = mask(var_v(1,:,:),((prob_u(0,:,:).lt.siglvl).or.(prob_v(0,:,:).lt.siglvl)),True)
copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
;copy_VarMeta(var(0,0,:,:),u_test(:,:))
;copy_VarMeta(var(0,0,:,:),v_test(:,:))
printVarSummary(var)
printVarSummary(prob)
printVarSummary(u_test)

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
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    
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
    res@lbBoxLineThicknessF  = mp_thick-2

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
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/0.01/)
    ;res_sig@cnFillPattern        = "GMT_gray"
    ;res_sig@cnFillColors         = (/"grey30","grey50"/) ;
    
    res_sig@gsnCenterString      = ""
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "black"
    resvc@vcMinDistanceF             = 0.03
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    resvc@vcRefAnnoOrthogonalPosF    = -0.68           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "
    resvc@vcLineArrowHeadMinSizeF   = 0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.015

    resvc2          = resvc
    resvc2@vcLineArrowColor           = "green"
    resvc2@vcLineArrowThicknessF      = vc_thick
    resvc2@vcRefAnnoOn                = False
    resvc2@vcLineArrowHeadMinSizeF   = 0.02
    resvc2@vcLineArrowHeadMaxSizeF   = 0.02

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
    resp@gsnPaperOrientation = "landscape"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    ;resp@txFontHeightF = title_font 
    ;resp@txFontThicknessF = title_thick

    vcref         = (/3,0.5/)   ;the level of 850hpa    ctrl,nudg,diff
    minlvl = (/(/ 0 , -2.7 /),(/ -4.5, -0.9/)/)
    maxlvl = (/(/ 54,  2.7 /),(/  4.5,  0.9/)/)
    spacig = (/(/ 3 ,  0.3 /),(/  0.5,  0.1/)/)

do nc = 0, ncase-1,1
    do nv = 0, nvar-1,1
        nn = nv +1
        np = nvar*nc + nv
        res@gsnLeftString   = pre_case(nc) + nn + case(nc)
        res@gsnRightString  = draw_var(nv)   
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc,nv)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc,nv) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc,nv)       ; set min contour level
        plot(np)     = gsn_csm_contour_map_ce(wks,var(nc,nv,:,:),res)
        if(nc.gt.0) then
            plot2(nv)  = gsn_csm_contour(wks,prob(nc-1,nv,:,:),res_sig)
            overlay(plot(np),plot2(nv))
        end if
        
        if(nv.eq.0) then
            resvc@vcRefMagnitudeF       = vcref(nc)             ; define vector ref mag
            resvc@vcRefAnnoString2      = vcref(nc)+"kg/(s*m)"
            plot_uv(nc)  = gsn_csm_vector(wks,var_u(nc,:,:),var_v(nc,:,:),resvc)
            overlay(plot(np),plot_uv(nc))
        end if

        if(nv.eq.0.and.nc.gt.0) then
            resvc2@vcRefMagnitudeF       = vcref(nc)             ; define vector ref mag
            resvc2@vcRefAnnoString2      = vcref(nc)+"kg/(s*m)"
            plot_uv2 = gsn_csm_vector(wks,u_test,v_test,resvc2)
            overlay(plot(np),plot_uv2)
        end if

        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var(nc,nv,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nv,{lats:latn},{lonl:lonr})))
    end do 
end do
    ;resp@txString      = fig_title + months(nm) 
    gsn_panel(wks,plot,(/ncase,nvar/),resp)
end

