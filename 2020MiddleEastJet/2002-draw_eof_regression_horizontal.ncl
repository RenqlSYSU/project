;***************************************************************************
;Draw the regression pattern based on the first three eof_ts 
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;                                            by Ql_Ren
;                                           2020/02/25
;******************************************************************************
begin
draw_var = (/"uwnd"/)  ;this is the variable used to draw, can used as center string
var_name = (/"season_ave"/)  ;use to read variables
unit     = (/"m/s"/)  ;this can be used as right string
nv = 0
lev  = (/850,500,200,100/)
nl = 2

area = (/"ME","EA","NA"/);,"GB"
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 0  ,105 , 250/);, 0   ;30 ;20
lonr = (/ 70 ,175 , 320/);, 360 ;150;220
lon_spacing = 20
lat_spacing = 10
narea = dimsizes(area)
neval = 4
siglvl = 0.05
years = ispan(1979,2015,1)

path   = "/home/ys17-19/renql/project/2020MiddleEastJet"
filename   = path + "/data/eof_u_6kinds-" + area + ".nc"
fig_name   = path + "/fig/eof_regression_hor-" + draw_var(nv)
fig_out    = "pdf" ;"eps";"x11"
pre_case   = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i", "(j","(k","(l", "(m","(n","(o", "(p","(q","(r"/)
plot   = new((/18/),graphic)
line   = new((/18/),graphic)
sigl   = new((/18/),graphic)

;===============================================
;draw the figure
;=============================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font

do na = 0, narea-1,1
;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing 
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False

    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.2
    res@pmLabelBarParallelPosF       = 0.5
    ;res@pmLabelBarWidthF  = 0.7
    ;res@pmLabelBarHeightF = 0.1
    
;---------set the attribute of the contour-----------------
    res_line = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True ;False
    res_line@cnLineLabelFontHeightF = 0.015
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    res_line@gsnContourZeroLineThicknessF = 8
    levels = new((/10/),integer)
    levels(0) = 0
    levels(1:9) = ispan(5,85,10)
    res_line@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    res_line@cnLevels             = levels       ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

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
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait" ;
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Vertical" ;"Horizontal"
    resp@pmLabelBarWidthF  = 0.06
    resp@pmLabelBarHeightF = 0.8
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(nv.eq.0) then
    minlvl = (/-7 /)
    maxlvl = (/ 7 /)
    spacig = (/ 1 /)
end if

f  = addfile(filename(na),"r")
ts = f->eof_ts
eof0 := f->eof0
eof1 := f->eof1
eof2 := f->eof2
eof3 := f->eof3
eof4 := f->eof4
eof5 := f->eof5
do ne = 0, neval-1, 1
    ts(0,ne,:) = ts(0,ne,:)/sqrt(eof0@eval(ne))
    ts(1,ne,:) = ts(1,ne,:)/sqrt(eof1@eval(ne))
    ts(2,ne,:) = ts(2,ne,:)/sqrt(eof2@eval(ne))
    ts(3,ne,:) = ts(3,ne,:)/sqrt(eof3@eval(ne))
    ts(4,ne,:) = ts(4,ne,:)/sqrt(eof4@eval(ne))
    ts(5,ne,:) = ts(5,ne,:)/sqrt(eof5@eval(ne))
end do

;f2 = addfile(path + "/data/ERA-Interim_DJF_" + draw_var(nv) + "-year.nc","r")
f2 = addfile(path + "/data/ERA-Interim_DJF_uwnd-year.nc","r")
vars := f2->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
u := dim_avg_n_Wrap(vars,0)
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)

var := regCoef_n(ts,vars,2,0)
copy_VarCoords(vars(0,:,:), var(0,0,:,:))
printVarSummary(var)
tval  := onedtond(var@tval ,(/6,neval,nlat,nlon/))
df    := onedtond(var@nptxy,(/6,neval,nlat,nlon/))
b = tval 
b = 0.5
prob  := betainc( df/(df+tval^2), df/2.0, b )
prob  := mask(prob,prob.lt.siglvl,True)
copy_VarCoords(vars(0,:,:), prob(0,0,:,:))
print(max(tval))
print(min(tval))
print(max(df))
print(min(df))

res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nk = 0, 5, 1
do ne = 0, neval-2, 1
    nn  = nk + 1
    nev = ne + 1
    np = nk*(neval-1) + ne 
    res@gsnLeftString   = pre_case(np) + ") " + nn + " EOF" + nev
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")" 
    if(nv.gt.10) then
    res@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
    else
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig        ; set contour interval
    res@cnMaxLevelValF        := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl       ; set min contour level
    end if
    plot(np)    = gsn_csm_contour_map_ce(wks,var(nk,ne,:,:),res)
    sigl(np)    = gsn_csm_contour(wks,prob(nk,ne,:,:),res_sig)
    overlay(plot(np),sigl(np))
    line(np)    = gsn_csm_contour(wks,u,res_line)
    overlay(plot(np),line(np))
    print(max(var(nk,ne,:,:)))
    print(min(var(nk,ne,:,:)))
end do 
end do 
    resp@txString = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"  
    gsn_panel(wks,plot,(/6,3/),resp)
end do ;area 
end

