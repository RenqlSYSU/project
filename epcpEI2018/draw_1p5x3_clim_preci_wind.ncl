;*********************************************************************
;filename:      draw_1p5x2_clim_preci_ano.ncl
;1 calculate the climatology of the var in a cycle of two years
;2 calc the difference of pSSTclim and pSSTpCP to get the preci anomaly
;3 draw the plot
;                                                     by Ql Ren
;                                                    2018/07/04
;*********************************************************************************
begin
;the path of the file
path    = "/HOME/sysu_hjkx_ys/WORKSPACE/renql/F/data"
case    = (/"pSSTclim","pSSTpEP","pSSTpCP"/)
file_in  = path + "/" + case + "/" + case + ".cam.h0.YEAR."
fig_name   = path + "/fig/" + (/"1p5x3_clim_preci_wnd"/)
fig_out    ="pdf";"x11"
unit       = "mm/day"   ;this can be used as right string 
draw_var = (/"preci"/)  ;this is the variable used to draw, can used as center string
var_name = (/"U","V","PRECC","PRECL"/) ;this is the var_name in data file
season   = (/"(a)JJA(0)","(b)SON(0)","(c)DJF(1)","(d)MAM(1)","(e)JJA(1)"/)

;area: east Asia
lats = 0 ;-20
latn = 50  ;60
lonl = 70  ;20
lonr = 180 ;220
lat_spac = 20
lon_spac = 40
lev  = (/850,500,200/)  ;can used as right string
nl   = 0
siglvl = 0.05
years  = ispan(1979,2004,1)

f    = addfile(file_in(0) + years(0) + ".monthly." + var_name(0) + ".nc","r")
vars = f->$var_name(0)$(0,{lev(nl)},{lats:latn},{lonl:lonr})
ncase  = dimsizes(case)
nsample = 13
nmonth  = 2*12
nseason = dimsizes(season)
nvar    = dimsizes(var_name)
nlat    = dimsizes(vars&lat) 
nlon    = dimsizes(vars&lon)
var     = new((/ncase,nvar-1,nsample,nmonth,nlat,nlon/),float)
var1    = new((/ncase,nvar-1,nsample,nseason,nlat,nlon/),float)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
    do nv = 0,nvar-3,1     ;read wind
    files := systemfunc("ls " + file_in(nc) + "*" + var_name(nv) + ".nc")
    do nt = 0,nsample-1,1
        f     := addfiles(files(2*nt:(2*nt+1)),"r")
        ListSetType(f,"cat")
        var(nc,nv,nt,:,:,:) = f[:]->$var_name(nv)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
    end do
    end do

    do nv = 2,nvar-1,1    ;read preci 
    files := systemfunc("ls " + file_in(nc) + "*" + var_name(nv) + ".nc")
    do nt = 0,nsample-1,1
        f     := addfiles(files(2*nt:(2*nt+1)),"r")
        ListSetType(f,"cat")
        var(nc,2,nt,:,:,:) = var(nc,2,nt,:,:,:) + f[:]->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    end do
    end do
end do
var(:,2,:,:,:,:) = var(:,2,:,:,:,:)*3600*24*1000 ;convert unit from m/s to mm/day

do ns = 0,nseason-1,1
    nms = ns*3 + 5
    nme = ns*3 + 7 
    var1(:,:,:,ns,:,:) = dim_avg_n(var(:,:,:,nms:nme,:,:),3) ;calc the season average
end do

var_ave  = dim_avg_n(var1,2)        ;get(ncase,nvar,nseason,nlat,nlon)
copy_VarMeta(vars,var_ave(0,0,0,:,:))
printVarSummary(var_ave)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;the number of the panel plot
    plot     = new((/nseason*ncase/),graphic)
    plotvc   = new((/nseason*ncase/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   vcref_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   vc_thick       = 3.0
   font_index     = 22

;draw the figure
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
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
    res@gsnMajorLatSpacing = lat_spac
    res@gsnMajorLonSpacing = lon_spac
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
    
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
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
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 1         ; set min contour level
    res@cnMaxLevelValF        = 16        ; set max contour level
    res@cnLevelSpacingF       = 1      ; set contour interval

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

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
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/0.01/)
    ;res_sig@cnFillPattern        = "GMT_gray"
    ;res_sig@cnFillColors         = (/"grey30","grey50"/) ;
    
    res_sig@gsnCenterString      = ""
    res_sig@gsnRightStringFontHeightF = subtitle_font
    res_sig@gsnLeftStringFontHeightF  = subtitle_font
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
    resvc@vcMinDistanceF             = 0.035
    resvc@vcRefMagnitudeF            = 10.0             ; define vector ref mag
    resvc@vcRefAnnoString2           = "10 m/s"
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = -0.67           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "
    resvc@vcLineArrowHeadMinSizeF   = 0.01
    resvc@vcLineArrowHeadMaxSizeF   = 0.01

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString       = " "
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString  = draw_var + " siglvl: " + siglvl

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

do ns = 0,nseason-1,1
do nc = 0,ncase-1,1
    np = ncase*ns+nc
    res@gsnLeftString   = season(ns) + " " + case(nc)
    plot(np)     = gsn_csm_contour_map_ce(wks,var_ave(nc,2,ns,:,:),res)
    plotvc(np)   = gsn_csm_vector(wks,var_ave(nc,0,ns,:,:),var_ave(nc,1,ns,:,:),resvc)
    overlay(plot(np),plotvc(np))
    print(max(var_ave(nc,2,ns,{lats:latn},{lonl:lonr})))
    print(min(var_ave(nc,2,ns,{lats:latn},{lonl:lonr})))
end do 
end do 

    resp@txString  = "climate" 
    gsn_panel(wks,plot,(/5,3/),resp)
end

