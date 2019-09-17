;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation(or EKE) in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;
;Modification
;2019/07/12: 1 add line of 200hPa u=20m/s and axis of westerly jet
;            2 add contour of 500hPa Z3 (geopotential height)
;            3 add case of noTP
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
draw_var = (/"intEKE"     ,"U"  ,"Z3"/)  ;this is the variable used to draw, can used as center string
unit     = (/"10~S~4~N~ J","m/s","m" /)  ;this can be used as right string 
;var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;months   = (/"DJF","MAM","JJA","SON"/)
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
lev = (/850,500,200/)
u_nl = 2 ;used to read westerly jet
z_nl = 1 ;used to read geopotential height
nv = 0
imonth = ispan(0,11,1)

;the path of the file
path    = "/home/ys17-19/renql/project/TP_NUDG"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.1) then
    case1   = (/"TP_CTRL","TP_CR"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
else
    case1   = (/"CTRL","NUDG6h"/);,"NUDG24h"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
end if

filename  = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" + draw_var(nv) + ".nc"
ufilename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" + draw_var(1) + ".nc"
zfilename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" + draw_var(2) + ".nc"
fig_name = path + "/analysis/fig/1p3x4_" + var_name(0) + "_" + draw_var(1) + "_" + draw_var(2)
;fig_title  = "Climate " + draw_var(nv) + " (" + unit(nv) + ")"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = 5 ;20 ;-20
latn = 55  ;55 ;60
lonl = 30  ;30 ;20
lonr = 150 ;150;220
siglvl = 0.05
lat_spacing = 15
lon_spacing = 30

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case1)
nmonth = dimsizes(imonth)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
u    = new((/ncase  ,nmonth,nlat,nlon/),float)
z    = new((/ncase  ,nmonth,nlat,nlon/),float)
var  = new((/ncase  ,nmonth,nlat,nlon/),float)
vari = new((/ncase  ,nmonth,nlat,nlon/),float)
prob = new((/ncase-1,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/nmonth/),graphic)
    gtop     = new((/nmonth/),graphic)
    line     = new((/nmonth/),graphic)
    plotu1   = new((/nmonth/),graphic)
    plotz1   = new((/nmonth/),graphic)
    axisu1   = new((/nmonth/),graphic)
    axisz1   = new((/nmonth/),graphic)
    plotu2   = new((/nmonth/),graphic)
    plotz2   = new((/nmonth/),graphic)
    axisu2   = new((/nmonth/),graphic)
    axisz2   = new((/nmonth/),graphic)
    dash_pattern = 16

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
   mp_thick       = 1.5
   font_index     = 22

;================================================================
;read data
;====================================================================
axlon = phis&lon
axlat = new((/ncase,nmonth,nlon/),double)
lat1 = phis&lat
print(lat1)
do nf = 0 ,ncase-1, 1
    z_f = addfile(zfilename(nf),"r")
    z(nf,:,:,:)  = z_f->$var_name(0)$(imonth,{lev(z_nl)},{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    
    u_f = addfile(ufilename(nf),"r")
    u(nf,:,:,:)  = u_f->$var_name(0)$(imonth,{lev(u_nl)},{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    do nm = 0, nmonth-1,1
    do inln = 0, nlon-1,1
        lat_ind = maxind(u(nf,nm,{5:50},inln)); + ind(lat1.eq.20)
        axlat(nf,nm,inln) = lat1(lat_ind)
    end do
    end do
    printVarSummary(axlat)
    
    f = addfile(filename(nf),"r")
    var(nf,:,:,:)  = f->$var_name(0)$(imonth,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,:,:,:) = f->$var_name(1)$(imonth,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
end do

do nf = 0 ,ncase-2,1
prob(nf,:,:,:) = ttest(var(0,:,:,:),vari(0,:,:,:),nsample,var(nf+1,:,:,:),vari(nf+1,:,:,:),nsample,True,False)
var(nf+1,:,:,:)= var(nf+1,:,:,:) - var(0,:,:,:)   ;get the difference between CTRL and NUDG
end do

var = var/10000

;var(1:2,:,:,:) = mask(var(1:2,:,:,:),prob.lt.siglvl,True)
var(1,:,:,:) = mask(var(1,:,:,:),prob(0,:,:,:).lt.siglvl,True)
prob = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
printVarSummary(var)
printVarSummary(u)
printVarSummary(z)
printVarSummary(prob)

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
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
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
    shadres = res
    shadres@cnFillOn         = True
    shadres@cnLinesOn        = False ;True
    shadres@cnLineLabelsOn   = False
    shadres@cnInfoLabelOn    = False
    ;shadres@cnLineThicknessF = 0
    ;shadres@gsnContourZeroLineThicknessF = 4.0
    ;shadres@cnLineColor      = "blue"
   
    ;shadres@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;shadres@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;shadres@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;shadres@cnLevels             = ispan(0,60,10)
    ;shadres@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;shadres@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;shadres@cnMinLevelValF        = 1          ; set min contour level
    ;shadres@cnMaxLevelValF        = 13 ;32        ; set max contour level
    ;shadres@cnLevelSpacingF       = 1        ; set contour interval

    shadres@gsnSpreadColors      = False  ;True              ; turn of the full color map
    shadres@gsnSpreadColorStart  = 11                ; start at color 17
    shadres@gsnSpreadColorEnd    = 2                ; end at color 200
    
    shadres@lbLabelBarOn         = True   ; the label bar
    shadres@lbOrientation        = "Horizontal"
    shadres@lbLabelFontHeightF   = tick_font              ; make labels smaller
    shadres@lbLabelFont          = font_index              ; make labels smaller
    shadres@lbLabelStride        = 1
    shadres@lbBoxLineThicknessF  = mp_thick-1

    shadres@pmLabelBarZone               = 0
    shadres@pmLabelBarOrthogonalPosF     = 0.15
    shadres@pmLabelBarParallelPosF       = 0.5

;----------Set res line-----------
    uline                  = True
    uline@gsnDraw          = False
    uline@gsnFrame         = False
    uline@gsnAddCyclic     = False
    uline@cnFillOn         = False       ; Fill contours
    uline@cnLinesOn        = True      ; turn off contour lines
    uline@cnLineLabelsOn   = False      ; turn off contour line labels
    uline@cnInfoLabelOn    = False
    uline@cnLineDashPattern = 0
    uline@cnLineThicknessF = 2.5
    uline@cnMonoLineColor  = True ;False
    uline@cnLineColor      = "blue"
    uline@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    uline@cnLevels         = (/20/)
    uline@lbLabelBarOn     = False   ; the label bar
    uline@gsnCenterString  = ""
    uline@gsnLeftString    = " "
    uline@gsnRightString   = " "

    uline2 = uline
    uline2@cnLineDashPattern = dash_pattern

    uaxis                      = True
    uaxis@gsLineColor          = "purple"
    uaxis@gsLineDashPattern    = 0
    uaxis@gsLineThicknessF     = 2.5

    uaxis2 = uaxis
    uaxis2@gsLineDashPattern    = dash_pattern

;----------Set res line for hgt-----------
    zline                  = True
    zline@gsnDraw          = False
    zline@gsnFrame         = False
    zline@gsnAddCyclic     = False
    zline@cnFillOn         = False       ; Fill contours
    zline@cnLinesOn        = True      ; turn off contour lines
    zline@cnLineLabelsOn   = False      ; turn off contour line labels
    zline@cnInfoLabelOn    = False
    zline@cnLineDashPattern = 0
    zline@cnLineThicknessF = 2.5
    zline@cnMonoLineColor  = True ;False
    zline@cnLineColor      = "orange" 
    zline@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    zline@cnLevels         = (/5860/)
    ;zline@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;zline@cnMinLevelValF        = 5000  ;4          ; set min contour level
    ;zline@cnMaxLevelValF        = 5800 ;60 ;32        ; set max contour level
    ;zline@cnLevelSpacingF       = 50 ;        ; set contour interval
    zline@gsnCenterString  = ""
    zline@gsnLeftString    = " "
    zline@gsnRightString   = " "

    zline2 = zline
    zline2@cnLineDashPattern = dash_pattern

    zaxis = zline
    zaxis@cnLineDashPattern = 0
    zaxis@cnLineThicknessF = 2.5
    zaxis@cnLineColor      = "red" 
    zaxis@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    zaxis@cnLevels         = (/0/)

    zaxis2 = zaxis
    zaxis2@cnLineDashPattern = dash_pattern

;--------set the attribute of polygon to draw the nudging  box--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 2.0
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

if(nv.eq.0) then ;EKE
    minlvl = (/ 2 ,-16,-16/);-8,-8,-8,-8
    maxlvl = (/50 , 16, 16/);8 ,8 ,8 ,8 
    spacig = (/ 3 , 2 , 2 /);1 ,1 ,1 ,1 
end if

    nc = 1
    res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
    do nm = 0,nmonth-1,1
        res@gsnLeftString   = pre_case(nm) + ") " + months(nm)
        res@gsnRightString  = "200U & 500hgt" 
        res@gsnSpreadColors       = False              ; turn of the full color map
        if(nv.gt.10) then
        shadres@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
        else
        shadres@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        shadres@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        shadres@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        shadres@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        end if
        ;plot(nm)    = gsn_csm_contour_map_ce(wks,var(nc,nm,:,:),res)
        plot(nm) = gsn_csm_map_ce(wks,res)
        gtop(nm) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(nm),gtop(nm))
        line(nm) = gsn_add_polyline( wks, plot(nm), lonx, latx, lnres )
        plotu1(nm)   = gsn_csm_contour(wks,u(0,nm,:,:),uline)
        plotz1(nm)   = gsn_csm_contour(wks,z(0,nm,:,:),zline)
        plotu2(nm)   = gsn_csm_contour(wks,u(1,nm,:,:),uline2)
        plotz2(nm)   = gsn_csm_contour(wks,z(1,nm,:,:),zline2)
        overlay(plot(nm),plotu1(nm))
        overlay(plot(nm),plotz1(nm))
        overlay(plot(nm),plotu2(nm))
        overlay(plot(nm),plotz2(nm))
        axisz1(nm)   = gsn_csm_contour(wks,u(0,nm,:,:),zaxis)
        axisz2(nm)   = gsn_csm_contour(wks,u(1,nm,:,:),zaxis2)
        overlay(plot(nm),axisz1(nm))
        overlay(plot(nm),axisz2(nm))
        axisu1(nm) = gsn_add_polyline( wks, plot(nm), axlon, axlat(0,nm,:), uaxis )
        axisu2(nm) = gsn_add_polyline( wks, plot(nm), axlon, axlat(1,nm,:), uaxis2 )
        print(max(var(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nm,{lats:latn},{lonl:lonr})))
        print(max(u(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(u(nc,nm,{lats:latn},{lonl:lonr})))
        print(max(z(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(z(nc,nm,{lats:latn},{lonl:lonr})))
    end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/3,4/),resp)
end

