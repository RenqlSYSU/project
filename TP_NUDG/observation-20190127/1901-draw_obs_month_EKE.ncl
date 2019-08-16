;*********************************************************************
;draw shaded monthly EKE or other variables
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
draw_var = (/"preci" ,"EKE"        ,"U"  ,"500hgt"/)  ;this is the variable used to draw, can used as center string
unit     = (/"mm/day","10~S~4~N~ J","m/s","m~S~2~N~/s~S~2~N~"/)    ;this can be used as right string 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
;var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
lev = (/850,500,200/)
nl = 1
nv = 1

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5 

;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127"
path1  = "/home/ys17-19/data-observation/ERA-interim/"
filename = path + "/mdata/ERA-Interim_Clim_month_ave_int" + draw_var(nv) + ".nc"
fig_name = path + "/fig/ERA-Interim_" + draw_var(nv) + "_" + var_name(0)
;filename = path + "/mdata/ERA-Interim_Clim_month_ave_int" + draw_var(nv) + "-" + cb + "-" + ca + "day.nc" 
;fig_name = path + "/fig/ERA-Interim_" + draw_var(nv) + "_" + var_name(0) + "-" + cb + "-" + ca + "day"
gtopfile = path1 + "gtopo_1.5x1.5.nc"
fig_title  = "Climate " + draw_var(nv) + " (" + unit(nv) + ")"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats =  -90 ;-15 ;20 ;-20
latn =  90 ;55  ;55 ;60
lonl =  0  ;30  ;30 ;20
lonr = 360 ;150 ;150;220
lat_spacing = 45
lon_spacing = 60

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats:latn},{lonl:lonr}))
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nmonth = dimsizes(months)
nlat   = dimsizes(phis&latitude)
nlon   = dimsizes(phis&longitude)
var  = new((/nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/nmonth/),graphic)
    gtop     = new((/nmonth/),graphic)
    line     = new((/nmonth/),graphic)

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
f = addfile(filename,"r")
var  = f->$var_name(0)$(:,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)

if(nv.eq.1) then
    var = var/10000
end if

printVarSummary(var)

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
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = 3.0
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

;---------set the attribute of the contour-----------------
    res_line = res
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = False
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.0
;    res_line@gsnContourZeroLineThicknessF = 3.0
    res_line@cnLineColor      = "blue"
   
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res_line@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res_line@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res_line@cnLevels             = ispan(0,60,10)
    ;res_line@cnFillColors         = (/-1,20,47,57,94,127,152/)
    ;res_line@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res_line@cnMinLevelValF        = 1          ; set min contour level
    ;res_line@cnMaxLevelValF        = 13 ;32        ; set max contour level
    res_line@cnLevelSpacingF       = 10        ; set contour interval

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(nv.eq.0) then ;preci
    minlvl = (/ 1 /);-8,-8,-8,-8
    maxlvl = (/16 /);8 ,8 ,8 ,8 
    spacig = (/ 1 /);1 ,1 ,1 ,1 
end if

if(nv.eq.1) then ;EKE
    minlvl = (/ 2 /);-8,-8,-8,-8
    maxlvl = (/98 /);8 ,8 ,8 ,8 
    spacig = (/6/);1 ,1 ,1 ,1 
    ;maxlvl = (/82/);8 ,8 ,8 ,8 
    ;spacig = (/5/);1 ,1 ,1 ,1 
end if

    res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
do nm = 0,nmonth-1,1
    res@gsnLeftString   = pre_case(nm) + ") Obs " + months(nm)
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")"
    res_line@gsnLeftString   = pre_case(nm) + ") Obs " + months(nm)
    res_line@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")"
    res@gsnSpreadColors       = False              ; turn of the full color map
    if(nv.gt.2) then
    res@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
    else
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig        ; set contour interval
    res@cnMaxLevelValF        := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl       ; set min contour level
    end if
    plot(nm) = gsn_csm_contour_map_ce(wks,var(nm,:,:),res) ;draw shaded figure
    ;plot(nm) = gsn_csm_contour_map_ce(wks,var(nm,:,:),res_line) ;draw contour figure
    gtop(nm) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(nm),gtop(nm))
    print(max(var(nm,{lats:latn},{lonl:lonr})))
    print(min(var(nm,{lats:latn},{lonl:lonr})))
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/2,2/),resp)
end

