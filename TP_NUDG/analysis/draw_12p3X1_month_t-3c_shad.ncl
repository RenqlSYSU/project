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
draw_var = (/"T","OMEGA","Q","Z3","DTCOND"/)  ;this is the variable used to draw, can used as center string
unit     = (/"~S~o~N~C","Pa/s","g/kg","km","K/day"/) ;"K"   ;this can be used as right string 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months   = (/"DJF","MAM","JJA","SON"/)
lev  = (/850,500,200/)  ;can used as right string
nl = 2
nv = 4

;the path of the file
;path     = "/users/yangsong3/renql/project/TP_NUDG"
path     = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_ave_" + draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/4p3X1_" + draw_var(nv) + "_" + var_name(0) + "-3c-" + lev(nl)
fig_title  = "Climate Seasonal " + draw_var(nv) + " in " ;+season
    ;case   = (/"AMIP_CTRL","(CTRL-NUDG24h)/CTRL%","(CTRL-NUDG6h)/CTRL%"/)
    case   = (/"AMIP_CTRL","NUDG24h-CTRL","NUDG6h-CTRL"/)
fig_out    ="pdf";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
siglvl = 0.05

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
nmonth = dimsizes(months)
nlat   = getfilevardimsizes(f,"lat")
nlon   = getfilevardimsizes(f,"lon")
nlev   = getfilevardimsizes(f,"lev")
var  = new((/ncase  ,nmonth,nlat,nlon/),float)
vari = new((/ncase  ,nmonth,nlat,nlon/),float)
prob = new((/ncase-1,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/3/),graphic)
    plot2    = new((/2/),graphic)
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
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf),"r")
    var(nf,:,:,:)  = f->$var_name(0)$(:,{lev(nl)},:,:) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,:,:,:) = f->$var_name(1)$(:,{lev(nl)},:,:);get an array (ncase,nmonth,nlat,nlon)
end do

if(nv.eq.0) then
    var      = var - 273.15 ;convert unit from K to c
    var_temp = where(var(0,:,:,:).eq.0,var@_FillValue,var(0,:,:,:))
    do nf = 0 ,ncase-2,1
        prob(nf,:,:,:) = ttest(var(0,:,:,:)+273.15,vari(0,:,:,:),27,var(nf+1,:,:,:)+273.15,vari(nf+1,:,:,:),27,True,False)
        var(nf+1,:,:,:)= (var(0,:,:,:) - var(nf+1,:,:,:))/abs(var_temp)*100 ;get the difference between CTRL and NUDG
    end do
end if

if(nv.eq.1) then
    var_temp = where(var(0,:,:,:).eq.0,var@_FillValue,var(0,:,:,:))
    do nf = 0 ,ncase-2,1
        prob(nf,:,:,:) = ttest(var(0,:,:,:),vari(0,:,:,:),27,var(nf+1,:,:,:),vari(nf+1,:,:,:),27,True,False)
        var(nf+1,:,:,:)= (var(0,:,:,:) - var(nf+1,:,:,:))/abs(var_temp)*100 ;get the difference between CTRL and NUDG
    end do
end if

if(nv.eq.2) then
    var      = var*1000      ;convert unit from kg/kg to g/kg
    var_temp = where(var(0,:,:,:).eq.0,var@_FillValue,var(0,:,:,:))
    do nf = 0 ,ncase-2,1
        prob(nf,:,:,:) = ttest(var(0,:,:,:)/1000,vari(0,:,:,:),27,var(nf+1,:,:,:)/1000,vari(nf+1,:,:,:),27,True,False)
        var(nf+1,:,:,:)= (var(0,:,:,:) - var(nf+1,:,:,:))/abs(var_temp)*100 ;get the difference between CTRL and NUDG
    end do
end if

if(nv.eq.3) then
;    var      = var/1000      ;convert unit from m to km
    var_temp = where(var(0,:,:,:).eq.0,var@_FillValue,var(0,:,:,:))
    do nf = 0 ,ncase-2,1
        prob(nf,:,:,:) = ttest(var(0,:,:,:)*1000,vari(0,:,:,:),27,var(nf+1,:,:,:)*1000,vari(nf+1,:,:,:),27,True,False)
        var(nf+1,:,:,:)= (var(0,:,:,:) - var(nf+1,:,:,:))/abs(var_temp)*100 ;get the difference between CTRL and NUDG
    end do
end if

if(nv.eq.4) then
    var = var*3600*24 ;convert unit from K/s to K/day
    do nf = 0 ,ncase-2,1
        prob(nf,:,:,:) = ttest(var(0,:,:,:)/3600/24,vari(0,:,:,:),27,var(nf+1,:,:,:)/3600/24,vari(nf+1,:,:,:),27,True,False)
        var(nf+1,:,:,:)= var(nf+1,:,:,:) - var(0,:,:,:)
    end do
end if

prob = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
printVarSummary(var)
printVarSummary(prob)
print(var(1,1,3,6))

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
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1

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
    res_sig@cnFillScaleF     = 1.
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

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 4.0
    ;lnres@gsLineLabelString    = "TP"

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

if(nv.eq.0) then
    minlvl = (/(/10  ,-9  ,-9  /),(/-18 ,-9  ,-9  /),(/ -64,-4.5,-4.5/)/)  ;850,500,200hpa
    maxlvl = (/(/28  , 9  , 9  /),(/  0 , 9  , 9  /),(/ -48, 4.5, 4.5/)/)
    spacig = (/(/ 1  , 1  , 1  /),(/  1 , 1  , 1  /),(/   1, 0.5, 0.5/)/) 
end if

if(nv.eq.1) then
    minlvl = (/(/-0.09,-90 ,-90 /),(/-0.09,-90 ,-90 /),(/-0.09,-90 ,-90 /)/)  ;850,500,200hpa
    maxlvl = (/(/ 0.09, 90 , 90 /),(/ 0.09, 90 , 90 /),(/ 0.09, 90 , 90 /)/)
    spacig = (/(/ 0.01, 10 , 10 /),(/ 0.01, 10 , 10 /),(/ 0.01, 10 , 10 /)/) 
end if

if(nv.eq.2) then
    minlvl = (/(/ 0  ,-18 ,-18 /),(/ 0  ,-18 ,-18 /),(/ 0  ,-18 ,-18 /)/)  ;850,500,200hpa
    maxlvl = (/(/18  ,18  , 18 /),(/4.5 ,18  , 18 /),(/0.18,18  , 18 /)/)
    spacig = (/(/ 1  , 2  ,  2 /),(/0.5 , 2  ,  2 /),(/0.01, 2  , 2  /)/) 
end if

if(nv.eq.3) then
    minlvl = (/(/1.25  ,-0.9 ,-0.9 /),(/5.0  ,-0.36,-0.36/),(/10.8,-0.36,-0.36/)/)  ;850,500,200hpa
    maxlvl = (/(/1.70  , 0.9 , 0.9 /),(/5.9  ,0.36 ,0.36 /),(/12.6,0.36 ,0.36 /)/)
    spacig = (/(/0.025 , 0.1 , 0.1 /),(/0.05 ,0.04 ,0.04 /),(/0.1 ,0.04 ,0.04 /)/) 
end if

if(nv.eq.4) then
    minlvl = (/(/ 0  ,-1.8 ,-1.8 /),(/ 0  ,-1.8 ,-1.8 /),(/ 0  ,-1.8 ,-1.8 /)/)
    maxlvl = (/(/ 6  ,1.8  , 1.8 /),(/ 6  ,1.8  , 1.8 /),(/ 6  ,1.8  , 1.8 /)/)
    spacig = (/(/0.5 ,0.2  , 0.2 /),(/0.5 ,0.2  , 0.2 /),(/0.5 ,0.2  , 0.2 /)/)
end if

do nm = 0 ,nmonth-1,1
    do nc = 0, ncase-1,1
        res@gsnLeftString   = case(nc)
        res@gsnSpreadColors       = False              ; turn of the full color map
        if(nv.eq.5) then
        res@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
        else
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       = spacig(nl,nc)        ; set contour interval
        res@cnMaxLevelValF        = maxlvl(nl,nc) ;32    ; set max contour level
        res@cnMinLevelValF        = minlvl(nl,nc)       ; set min contour level
        end if
        if(nc.gt.0) then
            res@gsnRightString  = "siglvl: " + siglvl
            plot(nc)     = gsn_csm_contour_map_ce(wks,var(nc,nm,:,:),res)
            plot2(nc-1)  = gsn_csm_contour(wks,prob(nc-1,nm,:,:),res_sig)
            overlay(plot(nc),plot2(nc-1))
        else
            res@gsnRightString   = unit(nv) + " " + lev(nl) + "hPa" ;draw_var(nv) 
            plot(nc)     = gsn_csm_contour_map_ce(wks,var(nc,nm,:,:),res)
        end if
        line(nc) = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
        print(max(var(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nm,{lats:latn},{lonl:lonr})))
    end do 

    resp@txString      = fig_title + months(nm) 
    gsn_panel(wks,plot,(/3,1/),resp)
end do ;for the month 
end

