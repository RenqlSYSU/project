;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
path    = "/home/ys17-19/renql/project/TP_NUDG"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG6h-CTRL"/) ;,"NG24h-CTRL"
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
end if
if(case_option.eq.1) then
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
end if
if(case_option.eq.2) then
    case1   = (/"F2000","NUDG"/)
    case2   = (/"F2000","NUDG-F2000"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 30
end if

;Q1 is calculated from temperature equation; q1 is calculated from (QRL,QRS,DTCOND,DTV)
draw_var = (/"intQ1"        ,"preci"       ,"OMEGA"           ,"Z3"             ,"intq1"     , \
             "T"            ,"500EKE"      ,"intEKE2-2.5-10day","intEKE2-10-90day","EKE"       , \
             "EKE-2.5-10day","EKE-10-90day","U"               ,"TS"             ,"PSL"       , \
             "precs" /) 
unit     = (/"W/kg"         ,"mm/day"      ,"10~s~-2~n~ pa/s" ,"gpm"            ,"W/kg"      , \
             "~S~o~N~C"     ,"10~S~4~N~J"  ,"10~S~4~N~J"      ,"10~S~4~N~J"     ,"J"         , \
             "J"            ,"J"           ,"m/s"             ,"~S~o~N~C"       ,"hPa"       , \
             "mm/day" /)   
draw_varname = (/"Q1"       ,"preci"       ,"omega"           ,"Z3"             ,"q1"        , \
                 "T"        ,"uppEKE"      ,"intEKEhgh"       ,"intEKElow"      ,"EKE"       , \
                 "EKEhgh"   ,"EKElow"      ,"U"               ,"TS"             ,"PSL"       , \
                 "snow" /) 
lev = (/850,500,200,100/)
nl = 1
nv = 15
diff = True ;False;
perc = False;True ;
vect_option = False;True ;
cont_option = False;True ;
jet_option  = 2 ;0 = no jet, 1 = jet axis, 2 = jet contour
sig_option  = True ;False;only work on shading 
cont_var  = (/"Z3" ,"intQ1","PSL" /)  ;this is the variable used to draw, can used as center string
cont_unit = (/"gpm","W/kg" ,"hPa"/)    ;this can be used as right string 
nvc = 2 
contsig  = True
barshare = False

var_option = 0 
if(var_option.eq.0) then 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
imonth = ispan(0,2,2)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
end if
if(var_option.eq.1) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;imonth = ispan(0,3,1);(/4,5,6,7/)
;imonth = ispan(4,7,1);(/4,5,6,7/)
;imonth = ispan(8,11,1);(/4,5,6,7/)
imonth = ispan(0,11,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
end if
if(var_option.eq.2) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"JF","MA","MJ","JA"/)
imonth = ispan(0,3,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_" ;+ draw_var(nv) + ".nc"
end if

fig_name = path + "/analysis/fig/" + case1(0) + "-2month_" + lev(nl) + "_" + draw_var(nv)
if(vect_option) then 
fig_name = path + "/analysis/fig/" + case1(0) + "-2month_" + lev(nl) + "wind_" + draw_var(nv)
end if
if(cont_option) then 
fig_name = path + "/analysis/fig/" + case1(0) + "-2month_" + lev(nl) + draw_var(nv) + "_" + cont_var(nvc)
end if

;pre_case   = (/"(e","(f","(g","(h","(e","(f","(g","(h","(i","(j","(k","(l"/)
pre_case   = (/"(a","(b","(d","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    = "pdf";"x11""eps";

;area: east Asia
lats = (/10 ,0  ,0  ,-70,-70,-70,10 ,-70, 0 /);-20
latn = (/89 ,70 ,70 , 0 , 0 , 0 ,55 ,70 ,60 /);60
lonl = (/0  ,150,270,30 ,150,270,90 ,0  ,40 /);20
lonr = (/360,270,390,150,270,390,150,360,150/);220
siglvl = 0.05
lat_spacing = 15
lon_spacing = 30
na = 0

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case1)
nmonth = dimsizes(imonth)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nvar = 4
var    = new((/ncase  ,nvar,nmonth,nlat,nlon/),float)
vari   = new((/ncase  ,nvar,nmonth,nlat,nlon/),float)
prob   = new((/ncase-1,nvar,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot   = new((/12/),graphic)
    plot2  = new((/12/),graphic)
    plot3  = new((/12/),graphic)
    wind   = new((/12/),graphic)
    wind2  = new((/12/),graphic)
    gtop   = new((/12/),graphic)
    line   = new((/12/),graphic)
    axis   = new((/12/),graphic)
    jetl   = new((/12/),graphic)
    sig    = new((/12/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;define the thickness and height of the title,tick and label
if(na.eq.0) then 
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   vcref_font     = 0.01
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 1.5
   font_index     = 22
else
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.0
   mp_thick       = 2.5
   font_index     = 22
end if

wks = gsn_open_wks(fig_out,fig_name)
do nl = 0,3,1
do nv = 10,11,1
;================================================================
;read data
;====================================================================
if(jet_option.ge.1) then 
axlon = phis&lon
axlat = new((/ncase,nmonth,nlon/),double)
lat1 = phis&lat
;print(lat1)
;print(ind(lat1.eq.20))
u    = new((/ncase  ,nmonth,nlat,nlon/),float)
do nf = 0 ,ncase-1, 1
    u_f = addfile(filename(nf)+"U.nc","r")
    u(nf,:,:,:)  = u_f->$var_name(0)$(imonth,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon):,:);
    do nm = 0, nmonth-1,1
    do inln = 0, nlon-1,1
        lat_ind = maxind(u(nf,nm,{10:48},inln))
        ;lat_ind = maxind(u(nf,nm,{20:50},inln)) + 37;ind(lat1.eq.20)
        axlat(nf,nm,inln) = lat1(lat_ind)
    end do
    end do
    printVarSummary(axlat)
end do
end if

do nf = 0 ,ncase-1, 1
if(vect_option) then 
inv = 0
f = addfile(filename(nf)+"U.nc","r")
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
inv = 1
f = addfile(filename(nf)+"V.nc","r")
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if

inv = 2
f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
if(nv.le.1.or.nv.eq.4.or.(nv.ge.6.and.nv.le.8).or.(nv.ge.13.and.nv.le.15)) then
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
else
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if

if(cont_option) then
inv = 3
f = addfile(filename(nf)+cont_var(nvc) + ".nc","r")
if(nvc.eq.1.or.nvc.eq.2) then
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
else
var(nf,inv,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
vari(nf,inv,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if
end if
end do

vari = where(abs(vari).le.0.000001,0.1,vari)
do nf = 0 ,ncase-2,1
prob(nf,:,:,:,:) = ttest(var(0,:,:,:,:),vari(0,:,:,:,:),nsample,var(nf+1,:,:,:,:),vari(nf+1,:,:,:,:),nsample,True,False)
end do

if(nvc.eq.2) then
    var(:,3,:,:,:) = var(:,3,:,:,:)/100 ;Pa -> hPa
end if

if(nv.eq.2) then
    var(:,2,:,:,:) = var(:,2,:,:,:)*100 ;omega*100
end if

if(nv.eq.4) then
    cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
    var(:,2,:,:,:) = var(:,2,:,:,:)*cp ;Q1,transform unit from K/S to W/kg
end if

if(nv.eq.5.or.nv.eq.13) then
    var(:,2,:,:,:) = var(:,2,:,:,:)-273.15 ;transform unit from K to c
end if

if(nv.ge.6.and.nv.le.8) then ;intEKE
    var(:,2,:,:,:) = var(:,2,:,:,:)/10000.0 ;transform unit from K to c
end if

if(nv.eq.14) then
    var(:,2,:,:,:) = var(:,2,:,:,:)/100 ;Pa -> hPa
end if

if(diff) then 
do nf = 0 ,ncase-2,1
var(nf+1,:,:,:,:)= var(nf+1,:,:,:,:) - var(0,:,:,:,:)   ;get the difference between CTRL and NUDG
if(perc) then 
var(0,2,:,:,:) = where(abs(var(0,2,:,:,:)).le.0.000001,0.1,var(0,2,:,:,:))
var(nf+1,2,:,:,:)= var(nf+1,2,:,:,:)*100.0/var(0,2,:,:,:)   ;get the difference between CTRL and NUDG
end if
end do

if(vect_option) then 
uvar = mask(var(1:(ncase-1),0,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),False)
vvar = mask(var(1:(ncase-1),1,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),False)
var(1:(ncase-1),0,:,:,:) = mask(var(1:(ncase-1),0,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),True)
var(1:(ncase-1),1,:,:,:) = mask(var(1:(ncase-1),1,:,:,:),((prob(:,0,:,:,:).lt.siglvl).or.(prob(:,1,:,:,:).lt.siglvl)),True)
copy_VarMeta(var(0,0,0,:,:),uvar(0,0,:,:))
copy_VarMeta(var(0,0,0,:,:),vvar(0,0,:,:))
end if
if(sig_option) then 
prob(:,2,:,:,:) = mask(prob(:,2,:,:,:),prob(:,2,:,:,:).lt.siglvl,True)
else
var(1:(ncase-1),2,:,:,:) = mask(var(1:(ncase-1),2,:,:,:),prob(:,2,:,:,:).lt.siglvl,True)
end if
if(cont_option.and.contsig) then 
var(1:(ncase-1),3,:,:,:) = mask(var(1:(ncase-1),3,:,:,:),prob(:,3,:,:,:).lt.siglvl,True)
end if
copy_VarMeta(var(0,0,0,:,:),prob(0,0,0,:,:))
printVarSummary(prob)
end if
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;    wks = gsn_open_wks(fig_out,fig_name)

;do na = 0, 5, 1
print("draw lats = " + lats(na) + ", latn = " + latn(na) + ", lonl = " + lonl(na) + ", lonr = " + lonr(na))
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "landscape" ;"portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
if(na.eq.2.or.na.eq.5) then
    res@mpCenterLonF     = 0.
else
    res@mpCenterLonF     = 180.
end if
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
    res@gsnAddCyclic     = False
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
   
    ;if(na.eq.0) then 
    ;res@lbLabelBarOn         = False   ; the label bar
    ;else
    ;res@lbLabelBarOn         = True   ; the label bar
    ;end if
    ;res@lbOrientation        = "Horizontal"
    ;res@lbLabelFontHeightF   = tick_font
    ;res@lbLabelFont          = font_index              ; make labels smaller
    ;res@lbLabelStride        = 1
    ;res@lbBoxLineThicknessF  = mp_thick-1
    ;res@pmLabelBarWidthF  = 0.65
    ;res@pmLabelBarHeightF = 0.08

    ;res@pmLabelBarZone               = 0
    ;res@pmLabelBarOrthogonalPosF     = 0.15
    ;res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
if(na.eq.0) then 
    resvc@vcMinDistanceF            = 0.0105
    resvc@vcRefLengthF              = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
else
    resvc@vcMinDistanceF            = 0.015 ;0.025
    resvc@vcRefLengthF              = 0.05 ;0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.01;0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.01;0.015
end if
    ;resvc@vcRefMagnitudeF            = 30.0             ; define vector ref mag
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.3           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    ;resvc@vcRefAnnoString2           = "30m/s"
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor           = "grey"

;----------Set res line for contour-----------
    res_line                  = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False       ; Fill contours
    res_line@cnLinesOn        = True      ; turn off contour lines
    res_line@cnLineLabelsOn   = False      ; turn off contour line labels
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnMonoLineColor  = True ;False
    if(jet_option.ge.1.or.vect_option) then
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    else
    res_line@cnLineColor      = "black"
    end if
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 2
    res_line@gsnContourPosLineDashPattern = 0
    res_line@cnLevelSelectionMode  = "AutomaticLevels"   ; set manual contour levels
    ;res_line@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res_line@cnLevelSpacingF       = 20 ;        ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;----------Set res line-----------
    res_jet                  = True
    res_jet@gsnDraw          = False
    res_jet@gsnFrame         = False
    res_jet@gsnAddCyclic     = False
    res_jet@cnFillOn         = False       ; Fill contours
    res_jet@cnLinesOn        = True      ; turn off contour lines
    res_jet@cnLineLabelsOn   = False      ; turn off contour line labels
    res_jet@cnInfoLabelOn    = False
    res_jet@cnLineThicknessF = 4
    res_jet@cnMonoLineColor  = True ;False
    if(vect_option) then
    res_jet@cnLineColor      = "blueviolet";"green4";"black"
    else
    res_jet@cnLineColor      = "black"
    end if
    if(jet_option.eq.1) then 
    res_jet@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_jet@cnLevels         = (/20/)
    end if
    if(jet_option.eq.2) then 
    u = mask(u, u.lt.10 , False) 
    res_jet@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_jet@cnLevelSpacingF      = 10        ; set contour interval
    end if
    res_jet@lbLabelBarOn         = False   ; the label bar
    res_jet@gsnCenterString      = ""
    res_jet@gsnLeftString    = " "
    res_jet@gsnRightString   = " "

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 2.0
    ;lnres@gsLineLabelString    = "TP"

    res_axis                      = True
    res_axis@gsLineColor          = "black"
    res_axis@gsLineDashPattern    = 1
    res_axis@gsLineThicknessF     = 3.5
    ;res_axis@gsLineLabelString    = "Westly"

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
    res_sig@cnFillDotSizeF   = 0.25
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

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
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = True
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.8 ;0.05
    resp@pmLabelBarHeightF = 0.05;0.8
    else
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Vertical";"Horizontal"
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.17
    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(na.eq.0) then 
    vcref         = (/(/10,5 ,5 /), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/10,5 ,5 /), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/30,10,10/), \  ;the level of 200hpa    ctrl,nudg,diff
                      (/20,10,10/)/)   ;the level of 300hpa    ctrl,nudg,diff
else
    vcref         = (/(/7 ,2 ,2 /), \  ;the level of 850hpa    ctrl,nudg,diff
                      (/7 ,2 ,2 /), \  ;the level of 700hpa    ctrl,nudg,diff
                      (/15,4 ,4 /), \  ;the level of 500hpa    ctrl,nudg,diff
                      (/15,4 ,4 /)/)   ;the level of 200hpa    ctrl,nudg,diff
end if

if(nvc.eq.0) then ;Z3
if(nl.eq.0) then ;850hPa
    cont_spacing = (/20,10/)
end if
if(nl.eq.1) then ;200hPa
    cont_spacing = (/80,10/)
end if
if(nl.eq.2) then ;200hPa
    cont_spacing = (/80,10/)
end if
end if

if(nvc.eq.1) then ;intQ1
    cont_spacing = (/20,1/)
end if

if(nvc.eq.2) then ;SLP
    cont_spacing = (/5,2/)
end if

if(nv.eq.0) 
    minlvl = (/0  ,-70,-140/)  ;q1 ;-80,-80
    maxlvl = (/300, 70, 140/)      ; 80, 80
    spacig = (/20 , 10,  20/)      ; 10, 10
end if
if(nv.eq.4) 
    minlvl = (/0  ,-160,-160/)  ;q1 ;-80,-80
    maxlvl = (/300, 160, 160/)      ; 80, 80
    spacig = (/20 ,  20,  20/)      ; 10, 10
end if

if(nv.eq.1) then ;preci
    minlvl = (/ 1 ,-1.4,-3.5/);-8,-8,-8,-8
    maxlvl = (/16 , 1.4, 3.5/);8 ,8 ,8 ,8 
    spacig = (/ 1 , 0.2, 0.5/);1 ,1 ,1 ,1 
end if

if(nv.eq.2) 
    minlvl = (/-8 ,-4 , -4/)  ;omega*100
    maxlvl = (/ 8 , 4 ,  4/)
    spacig = (/ 1 ,0.5,0.5/)
end if

if(nv.eq.3) ;Z3
if(nl.eq.0) then ;850 
    minlvl = (/1300,-35,-35/); -24,-24
    maxlvl = (/1525, 35, 35/);  24, 24
    spacig = (/ 15 , 5 , 5 /);   3,  3
end if
if(nl.eq.2) then  ;200
    minlvl = (/10900,-70,-70/) 
    maxlvl = (/12400, 70, 70/)
    spacig = (/  100, 10, 10/) 
end if
end if

if(nv.eq.5) ;T
if(nl.eq.0) then ;850 
    minlvl = (/-15,-3.5,-3.5/); -24,-24
    maxlvl = (/ 15, 3.5, 3.5/);  24, 24
    spacig = (/ 2 ,0.5 ,0.5 /);   3,  3
end if
if(nl.eq.1) then ;500 
    minlvl = (/-38,-3.5/); -24,-24
    maxlvl = (/-8 , 3.5/);  24, 24
    spacig = (/ 2 ,0.5 /);   3,  3
end if
if(nl.eq.2) then  ;200
    minlvl = (/-68,-3.5,-3.5/); -24,-24
    maxlvl = (/-50, 3.5, 3.5/);  24, 24
    spacig = (/1.2,0.5 ,0.5 /);   3,  3
end if
end if

if(nv.ge.6.and.nv.le.11) then ;intEKE
if(perc) then 
    minlvl = (/ 5 ,-70,-70/);-8,-8,-8,-8
    maxlvl = (/80 , 70, 70/);8 ,8 ,8 ,8 
    spacig = (/ 5 , 10, 10/);1 ,1 ,1 ,1 
else
    minlvl = (/10,-35/);-8,-8,-8,-8
    maxlvl = (/85, 35/);8 ,8 ,8 ,8 
    spacig = (/5 , 5 /);1 ,1 ,1 ,1 
end if
end if

if(nv.eq.12) 
    minlvl = (/ 0 ,-7 ,-7 /)  ;omega*100
    maxlvl = (/ 75, 7 , 7 /)
    spacig = (/ 5 , 1 , 1 /)
end if

if(nv.eq.13) 
    minlvl = (/-45 ,-3.5 ,-7 /)  ;omega*100
    maxlvl = (/ 30, 3.5 , 7 /)
    spacig = (/ 5 , 0.5 , 1 /)
end if

if(nv.eq.14) 
    minlvl = (/1000,-14,-7 /)  ;omega*100
    maxlvl = (/1030, 14, 7 /)
    spacig = (/ 2  , 2 , 1 /)
end if

if(nv.eq.15) then ;precs
    minlvl = (/0.5,-1.4,-3.5/);-8,-8,-8,-8
    maxlvl = (/ 8 , 1.4, 3.5/);8 ,8 ,8 ,8 
    spacig = (/0.5, 0.2, 0.5/);1 ,1 ,1 ,1 
end if

do nc = 0,1,1
;do nm = 0,nmonth-1,1
if((diff.and.nc.gt.0).or.nv.eq.2) then
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
   ; res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
   res@cnFillColors  := (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
else
    gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
   res@cnFillColors  := (/2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)
end if 
    do nm = 0,nmonth-1,1
        ;np = 2*nl+(nv-10) 
        np = nm
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(imonth(nm))
        if(nv.le.1.or.nv.eq.4.or.(nv.ge.6.and.nv.le.8).or.(nv.ge.13.and.nv.le.15)) then
        if(cont_option) then 
        res@gsnRightString  = draw_varname(nv) + " (" + unit(nv) + ") & " + cont_var(nvc) + " (" + cont_unit(nvc) + ") "
        else
        res@gsnRightString  = draw_varname(nv) + " (" + unit(nv) + ")" 
        end if
        else
        if(cont_option) then 
        res@gsnRightString  = lev(nl) + "hPa " + draw_varname(nv) + " (" + unit(nv) + ") & " + cont_var(nvc) + " (" + cont_unit(nvc) + ") "
        else
        res@gsnRightString  = lev(nl) + "hPa " + draw_varname(nv) + " (" + unit(nv) + ") "
        end if
        end if
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        plot(np)    = gsn_csm_contour_map_ce(wks,var(nc,2,nm,:,:),res)
        if(sig_option.and.nc.gt.0) then 
        sig(np) = gsn_csm_contour(wks,prob(nc-1,2,nm,:,:),res_sig)
        overlay(plot(np),sig(np))
        end if

        if(vect_option) then 
        resvc@vcRefMagnitudeF      = vcref(nl,nc)             ; define vector ref mag
        resvc@vcRefAnnoString2     = vcref(nl,nc)+"m/s"
        resvc2@vcRefMagnitudeF      = vcref(nl,nc)             ; define vector ref mag
        resvc2@vcRefAnnoString2     = vcref(nl,nc)+"m/s"
        if(diff.and.nc.gt.0)
        wind2(np)   = gsn_csm_vector(wks,uvar(nc-1,nm,:,:),vvar(nc-1,nm,:,:),resvc2)
        overlay(plot(np),wind2(np))
        end if
        wind(np)    = gsn_csm_vector(wks,var(nc,0,nm,:,:),var(nc,1,nm,:,:),resvc)
        overlay(plot(np),wind(np))
        end if
        
        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        
        if(cont_option)
        res_line@cnLevelSpacingF    = cont_spacing(nc) 
        ;res_line@cnLineColor      = "blueviolet";"green4";"black"
        plot2(np)   = gsn_csm_contour(wks,var(nc,3,nm,:,:),res_line)
        overlay(plot(np),plot2(np))
        ;if(nc.eq.1) then 
        ;res_line2   = res_line
        ;res_line2@cnLevelSpacingF  = cont_spacing(0)
        ;res_line2@cnLineColor      = "black"
        ;plot3(np)   = gsn_csm_contour(wks,var(0,3,nm,:,:),res_line2)
        ;overlay(plot(np),plot3(np))
        ;end if
        end if

        if(na.eq.0.or.na.eq.6) then
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        end if
        
        if(jet_option.ge.1) then 
        if(jet_option.eq.1) then 
        axis(np) = gsn_add_polyline( wks, plot(np), axlon, axlat(nc,nm,:), res_axis )
        end if
        jetl(np) = gsn_csm_contour(wks,u(nc,nm,:,:),res_jet)
        overlay(plot(np),jetl(np))
        end if
        print(max(var(nc,2,nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        print(min(var(nc,2,nm,{lats(na):latn(na)},{lonl(na):lonr(na)})))
        ;draw(plot(np))
        ;frame(wks)
    end do 
    ;gsn_panel(wks,plot,(/ncase,1/),resp)
;end do ;for the month 

    res1 = True
    res2 = True
    res1@gsnMaximize = True
    ;attach1 = gsn_attach_plots(plot(0),plot(1:3),res1,res2)
    ;attach2 = gsn_attach_plots(plot(4),plot(5:7),res1,res2)
    ;attach3 = gsn_attach_plots(plot(8),plot(9:11),res1,res2)
    ;gsn_panel(wks,plot(0:8:4),(/3,1/),resp)
    ;attach1 = gsn_attach_plots(plot(0),plot(1:2),res1,res2)
    ;attach2 = gsn_attach_plots(plot(3),plot(4:5),res1,res2)
    ;attach3 = gsn_attach_plots(plot(6),plot(7:8),res1,res2)
    ;attach4 = gsn_attach_plots(plot(9),plot(10:11),res1,res2)
    ;gsn_panel(wks,plot(0:9:3),(/4,1/),resp)
   ; attach1 = gsn_attach_plots(plot(0),plot(1),res1,res2)
   ; attach2 = gsn_attach_plots(plot(2),plot(3),res1,res2)
   ; attach3 = gsn_attach_plots(plot(4),plot(5),res1,res2)
   ; attach4 = gsn_attach_plots(plot(6),plot(7),res1,res2)
   ; attach5 = gsn_attach_plots(plot(8),plot(9),res1,res2)
   ; attach6 = gsn_attach_plots(plot(10),plot(11),res1,res2)
   ; gsn_panel(wks,plot(0:10:2),(/6,1/),resp)
;    resp@txString      = fig_title
    ;gsn_panel(wks,plot,(/2,2/),resp)
;    gsn_panel(wks,plot,(/2,1/),resp)
end do
end do
end do
    gsn_panel(wks,plot,(/4,2/),resp)
end

