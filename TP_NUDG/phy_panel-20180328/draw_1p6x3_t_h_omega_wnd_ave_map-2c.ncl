;*********************************************************************
;filename:      draw_1p2x3_wind_t_h_ave_map-2c.ncl
;input_data:    CTRL-Clim_month_ave_*.nc 
;
;1 read JJA data from four files season_ave(2,:,:) and season_vari(2,:,:)
;2 do the significance test
;3 define the levels then draw the plot 
;                                                     by Ql Ren
;                                                    2018/03/28
;*********************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.0
   mp_thick       = 4.0
   font_index     = 22

draw_var = (/"T(~S~o~N~C)","Z(gpm)","omega(10~S~-2~N~ Pa/s)"/)  ;this is the variable used to draw, can used as center string
draw_varname = (/"T","Z3","OMEGA","U","V"/)
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months   = (/"DJF","MAM","JJA","SON"/)
lev      = (/850,500,200/)
nm = 2

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
filename = path + "/analysis/mdata/" + (/"CTRL","NUDG6h"/) + "-Clim_month_ave_" ;+ draw_varname(nv) + ".nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/phy_panel-20180328/fig/1p6x3_t_h_omega_wnd_" + months(nm) + var_name(0) 
    case   = (/" CTRL"," NG6h-CTRL"/)
pre_case   = (/(/"(a","(b"/),(/"(c","(d"/),(/"(e","(f"/)/)
fig_out    = "eps";"pdf";"x11"
;fig_out    = "pdf";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case)
nlev   = dimsizes(lev)
nmonth = dimsizes(months)
nvar   = dimsizes(draw_var)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var    = new((/ncase  ,nvar,nlev,nlat,nlon/),float)
var_u  = new((/ncase       ,nlev,nlat,nlon/),float)
var_v  = new((/ncase       ,nlev,nlat,nlon/),float)
vari   = new((/ncase  ,nvar,nlev,nlat,nlon/),float)
vari_u = new((/ncase       ,nlev,nlat,nlon/),float)
vari_v = new((/ncase       ,nlev,nlat,nlon/),float)
prob   = new((/ncase-1,nvar,nlev,nlat,nlon/),float)
prob_u = new((/ncase-1     ,nlev,nlat,nlon/),float)
prob_v = new((/ncase-1     ,nlev,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase*nvar*nlev/),graphic)
    line     = new((/ncase*nvar*nlev/),graphic)
    gtop     = new((/ncase*nvar*nlev/),graphic)
    plot2    = new((/(ncase-1)*(nvar-1)*nlev/),graphic)
    plot_vc  = new((/2*nlev/),graphic)

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
do nv = 0 ,nvar-1,1
    f = addfile(filename(nf)+draw_varname(nv)+".nc","r")
     var(nf,nv,:,:,:)  = f->$var_name(0)$(nm,:,:,:) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,nv,:,:,:)  = f->$var_name(1)$(nm,:,:,:) ;get an array (ncase,nmonth,nlat,nlon)
end do
f = addfile(filename(nf)+draw_varname(3)+".nc","r")
 var_u(nf,:,:,:)  = f->$var_name(0)$(nm,:,:,:) ;get an array (ncase,nmonth,nlat,nlon)
vari_u(nf,:,:,:)  = f->$var_name(1)$(nm,:,:,:) ;get an array (ncase,nmonth,nlat,nlon)
f = addfile(filename(nf)+draw_varname(4)+".nc","r")
 var_v(nf,:,:,:)  = f->$var_name(0)$(nm,:,:,:) ;get an array (ncase,nmonth,nlat,nlon)
vari_v(nf,:,:,:)  = f->$var_name(1)$(nm,:,:,:) ;get an array (ncase,nmonth,nlat,nlon)
end do

do nf = 0 ,ncase-2,1
prob(nf,:,:,:,:) = ttest(var(0,:,:,:,:),vari(0,:,:,:,:),27,var(nf+1,:,:,:,:),vari(nf+1,:,:,:,:),27,True,False)
prob_u(nf,:,:,:) = ttest(var_u(0,:,:,:),vari_u(0,:,:,:),27,var_u(nf+1,:,:,:),vari_u(nf+1,:,:,:),27,True,False)
prob_v(nf,:,:,:) = ttest(var_v(0,:,:,:),vari_v(0,:,:,:),27,var_v(nf+1,:,:,:),vari_v(nf+1,:,:,:),27,True,False)
var(nf+1,:,:,:,:)= var(nf+1,:,:,:,:) - var(0,:,:,:,:) ;get the difference between CTRL and NUDG
var_u(nf+1,:,:,:)= var_u(nf+1,:,:,:) - var_u(0,:,:,:) ;get the difference between CTRL and NUDG
var_v(nf+1,:,:,:)= var_v(nf+1,:,:,:) - var_v(0,:,:,:) ;get the difference between CTRL and NUDG
end do

var(0,0,:,:,:) = var(0,0,:,:,:)-273.15      ;convert unit from K to C
var(:,2,:,:,:) = var(:,2,:,:,:)*100     ;convert unit from Pa/s to 10^-2 Pa/s

prob   = mask(prob,prob.lt.siglvl,True)
var(1,2,:,:,:) = mask(var(1,2,:,:,:),prob(0,2,:,:,:).lt.siglvl,True)
var_u(1,:,:,:) = mask(var_u(1,:,:,:),((prob_u(0,:,:,:).lt.siglvl).or.(prob_v(0,:,:,:).lt.siglvl)),True)
var_v(1,:,:,:) = mask(var_v(1,:,:,:),((prob_u(0,:,:,:).lt.siglvl).or.(prob_v(0,:,:,:).lt.siglvl)),True)
copy_VarMeta(var(0,0,:,:,:),prob(0,0,:,:,:))
printVarSummary(var)
printVarSummary(prob)

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
    
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    
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
    
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "grey0"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False   ; the label bar
    
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
    resvc@vcMinDistanceF             = 0.032
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.2           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@vcLineArrowHeadMinSizeF   = 0.012
    resvc@vcLineArrowHeadMaxSizeF   = 0.012
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

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

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait";"landscape"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    ;resp@txFontHeightF = title_font 
    ;resp@txFontThicknessF = title_thick

    vcref         = (/(/5 ,1/), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/5 ,2/), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/15,3/)/)   ;the level of 200hpa    ctrl,nudg,diff

do nl = 0,nlev-1,1
if(nl.eq.0) then
    minlvl = (/(/ 12,1410,-8  /),(/-0.8, -8 , -4.0/)/)
    maxlvl = (/(/ 27,1560, 8  /),(/ 0.8,  8 ,  4.0/)/)
    spacig = (/(/ 1 ,  10, 1  /),(/ 0.1,  1 ,  0.5/)/)
end if

if(nl.eq.1) then
    minlvl = (/(/-7.5,5730,-8  /),(/-0.8, -16, -4.0/)/)
    maxlvl = (/(/ 0  ,5880, 8  /),(/ 0.8,  16,  4.0/)/)
    spacig = (/(/0.5 ,  10, 1  /),(/ 0.1,  2 ,  0.5/)/)
end if

if(nl.eq.2) then
    minlvl = (/(/-55.5,12063,-8  /),(/-0.8, -24, -4.0/)/)
    maxlvl = (/(/-48  ,12558, 8  /),(/ 0.8,  24,  4.0/)/)
    spacig = (/(/0.5  ,  33 , 1  /),(/ 0.1,  3 ,  0.5/)/)
end if

do nc = 0, ncase-1,1
    do nv = 0, nvar-1,1
if(nc.eq.0.and.nv.lt.2) then
    res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
else 
    res@cnFillPalette = "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
end if 
        nn = nv +1
        np = nvar*ncase*nl + nvar*nc + nv
        res@gsnLeftString   = pre_case(nl,nc) + nn + ") " + case(nc)
        res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv)  
    if(nc.eq.1.and.nv.eq.2) then
        res@gsnRightString  = draw_var(nv)   
    end if
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc,nv)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc,nv) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc,nv)       ; set min contour level
        plot(np)     = gsn_csm_contour_map_ce(wks,var(nc,nv,nl,:,:),res)
        if(nv.eq.2) then
            resvc@vcRefMagnitudeF       = vcref(nl,nc)             ; define vector ref mag
            resvc@vcRefAnnoString2      = vcref(nl,nc)+"m/s"
            plot_vc(nc)  = gsn_csm_vector(wks,var_u(nc,nl,:,:),var_v(nc,nl,:,:),resvc)
            overlay(plot(np),plot_vc(nc))
        end if
        
        if(nc.gt.0.and.nv.lt.2) then
            plot2(nv)  = gsn_csm_contour(wks,prob(nc-1,nv,nl,:,:),res_sig)
            overlay(plot(np),plot2(nv))
        end if

        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var(nc,nv,nl,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nv,nl,{lats:latn},{lonl:lonr})))
    end do 
end do
end do
    gsn_panel(wks,plot,(/ncase*nlev,nvar/),resp)
end

