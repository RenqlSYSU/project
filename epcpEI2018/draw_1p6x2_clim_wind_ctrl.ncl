;*********************************************************************
;filename:      draw_1p6x2_clim_wind_ano.ncl
;1 calculate the climatology of the var in a years
;3 draw the plot
;                                                     by Ql Ren
;                                                    2018/07/04
;*********************************************************************************
begin
;area: east Asia
lats = -30 ;-20
latn = 50  ;60
lonl = 70  ;20
lonr = 280 ;220
lat_spac = 20
lon_spac = 40
lev  = (/850,200/)  ;can used as right string
nl = 0
vcref = (/10,30/)
siglvl = 0.05
years  = ispan(1979,2004,1)

;the path of the file
path    = "/HOME/sysu_hjkx_ys/WORKSPACE/renql/F/data"
case    = (/"pSSTclim"/)
file_in  = path + "/" + case + "/" + case + ".cam.h0"
fig_name   = path + "/fig/1p6x2_clim_wind_ctrl" + lev(nl) + "hPa-" + lonl + "-" + lonr
fig_out    ="pdf";"x11"
unit       = "m/s"   ;this can be used as right string 
draw_var = (/"wind"/)  ;this is the variable used to draw, can used as center string
var_name = (/"U","V"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

f    = addfile(file_in(0) + ".YEAR." + years(0) + ".monthly." + var_name(0) + ".nc","r")
vars = f->$var_name(0)$(0,{lev(nl)},{lats:latn},{lonl:lonr})
ncase  = dimsizes(case)
nsample = dimsizes(years)*2 
nmonth  = 12
nvar    = dimsizes(var_name)
nlat    = dimsizes(vars&lat) 
nlon    = dimsizes(vars&lon)
var     = new((/ncase,nvar,nsample,nmonth,nlat,nlon/),float)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
    do nv = 0,nvar-1,1     ;read wind
    files := systemfunc("ls " + file_in(nc) + "*" + var_name(nv) + ".nc")
    do nt = 0,nsample-1,1
        f     := addfile(files(nt),"r")
        var(nc,nv,nt,:,:,:) = f->$var_name(nv)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
    end do
    end do
end do

var_ave  = dim_avg_n(var,2)        ;get(ncase,nvar,nseason,nlat,nlon)
copy_VarMeta(vars,var_ave(0,0,0,:,:))
printVarSummary(var_ave)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;the number of the panel plot
    plot     = new((/nmonth/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   vcref_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   vc_thick       = 2.0
   font_index     = 22

;draw the figure
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"precip3_16lev")
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

;--------set the attribute of vector--------------
    res@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    res@vcVectorDrawOrder          = "PostDraw"
    res@vcLineArrowThicknessF      = vc_thick
    res@vcMonoLineArrowColor       = True
    res@vcLineArrowColor           = "black"
    res@vcMinDistanceF             = 0.035
    res@vcRefMagnitudeF            = vcref(nl)             ; define vector ref mag
    res@vcRefAnnoString2           = vcref(nl) + " m/s"
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
    res@vcLineArrowHeadMinSizeF   = 0.01
    res@vcLineArrowHeadMaxSizeF   = 0.01

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString       = " "
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString  = lev(nl) + "hPa " + draw_var + " (" + unit + ")"

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

do ns = 0,nmonth-1,1
    res@gsnLeftString   = months(ns) + " " + case 
    plot(ns)     = gsn_csm_vector_map_ce(wks,var_ave(0,0,ns,:,:),var_ave(0,1,ns,:,:),res)
    print(max(var_ave(0,0,ns,{lats:latn},{lonl:lonr})))
    print(min(var_ave(0,0,ns,{lats:latn},{lonl:lonr})))
end do

;    resp@txString  = "climate" 
    gsn_panel(wks,plot,(/6,2/),resp)
end

