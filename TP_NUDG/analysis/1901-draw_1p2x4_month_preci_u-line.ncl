;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;
;Modification:
;20190712: add the line of 200hPa u=20m/s and westerly axis
;
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
draw_var = (/"preci" ,"EKE"     ,"q1"  ,"U"  /)  ;this is the variable used to draw, can used as center string
unit     = (/"mm/day","10~S~4~N~ J","W/kg","m/s"/)  ;this can be used as right string 
;var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
;months   = (/"DJF","MAM","JJA","SON"/)
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
lev = (/850,500,200/)
nl = 2 ;used to read westerly jet
nv = 1
imonth = ispan(4,7,1);(/4,5,6,7/)
;imonth = ispan(0,3,1);(/4,5,6,7/)

;the path of the file
path    = "/home/ys17-19/renql/project/TP_NUDG"
;case1   = (/"TP_CTRL","TP_CR"/)
;case2   = (/"CTRL","noTP-CTRL"/)
;;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
;gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
;nsample = 10
case1   = (/"CTRL","NUDG6h","NUDG24h"/)
case2   = (/" CTRL"," NG6h-CTRL", " NG24h-CTRL"/)
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
nsample = 27

filename  = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" + draw_var(nv) + ".nc"
ufilename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
fig_name  = path + "/analysis/fig/1p2x4_" + draw_var(nv) + "_" + var_name(0) + "-2c" 
;fig_title  = "Seasonal " + draw_var(nv) + " variance in " ;+season
fig_title  = "Climate " + draw_var(nv) + " (" + unit(nv) + ")"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    = "eps";"pdf";"x11"

;area: east Asia
lats = -15 ;20 ;-20
latn = 55  ;55 ;60
lonl = 30  ;30 ;20
lonr = 150 ;150;220
siglvl = 0.05

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"North China","South China","TP_NUDG"/)
alats = (/30  ,20 ,25 /)
alatn = (/40  ,30 ,40 /)
alonl = (/107 ,107,70 /)
alonr = (/122 ,122,105/)

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case2)
nmonth = dimsizes(imonth)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
u    = new((/ncase  ,nmonth,nlat,nlon/),float)
var  = new((/ncase  ,nmonth,nlat,nlon/),float)
vari = new((/ncase  ,nmonth,nlat,nlon/),float)
prob = new((/ncase-1,nmonth,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase*nmonth/),graphic)
    plot2    = new((/ncase*nmonth/),graphic)
    gtop     = new((/ncase*nmonth/),graphic)
    line     = new((/ncase*nmonth/),graphic)
    axis     = new((/ncase*nmonth/),graphic)
    plot_sig = new((/(ncase-1)*nmonth/),graphic)

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
print(ind(lat1.eq.20))
do nf = 0 ,ncase-1, 1
    u_f = addfile(ufilename(nf),"r")
    u(nf,:,:,:)  = u_f->$var_name(0)$(imonth,{lev(nl)},{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    do nm = 0, nmonth-1,1
    do inln = 0, nlon-1,1
        lat_ind = maxind(u(nf,nm,{20:50},inln)) + 37;ind(lat1.eq.20)
        axlat(nf,nm,inln) = lat1(lat_ind)
    end do
    end do
    printVarSummary(axlat)
    
    f = addfile(filename(nf),"r")
    if(nv.eq.3) then
    fig_name  = path + "/analysis/fig/1p2x4_" + draw_var(nv) + "_" + lev(nl) + var_name(0) + "-2c" 
    fig_title = "Climate " + lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"
    var(nf,:,:,:)  = f->$var_name(0)$(imonth,{lev(nl)},{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,:,:,:) = f->$var_name(1)$(imonth,{lev(nl)},{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var(nf,:,:,:)  = f->$var_name(0)$(imonth,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,:,:,:) = f->$var_name(1)$(imonth,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end do

;------calc the area average q1 and print---------------
if(nv.eq.2) then
narea = dimsizes(area)                                                                              ;
do na = 0, narea-1,1                                                                                ;
    if(na.eq.2) then                                                                                ;
        var1 := var(:,:,{alats(na):alatn(na)},{alonl(na):alonr(na)})                                ;
        phis1:= conform(var1,phis({alats(na):alatn(na)},{alonl(na):alonr(na)}),(/2,3/))             ;
        var1 := mask(var1,phis1.gt.1500,True)                                                       ;
        printVarSummary(var1)                                                                       ;
        copy_VarMeta(var(:,:,{alats(na):alatn(na)},{alonl(na):alonr(na)}),var1)                     ;
        temp = wgt_areaave_Wrap(var1,1.0,1.0,0)                                                     ;
    else                                                                                            ;
        temp = wgt_areaave_Wrap(var(:,:,{alats(na):alatn(na)},{alonl(na):alonr(na)}),1.0,1.0,0)     ;
    end if                                                                                          ;
    do nm = 0, nmonth-1,1                                                                           ;
        print(months(imonth(nm)) + " " + area(na) + ": " + temp(:,nm))                              ;
    end do                                                                                          ;
end do          
end if
;----------------------------------------------                                                     
                                                                                                    
do nf = 0 ,ncase-2,1
prob(nf,:,:,:) = ttest(var(0,:,:,:),vari(0,:,:,:),nsample,var(nf+1,:,:,:),vari(nf+1,:,:,:),nsample,True,False)
var(nf+1,:,:,:)= var(nf+1,:,:,:) - var(0,:,:,:)   ;get the difference between CTRL and NUDG
end do

if(nv.eq.1) then ;EKE
    var = var/10000
end if

if(nv.eq.2) then   ;q1
    var = var*1000 ;transform unit from K/S to w/kg
end if

var(1:2,:,:,:) = mask(var(1:2,:,:,:),prob.lt.siglvl,True)      ;use for NUDG
;var(1,:,:,:) = mask(var(1,:,:,:),prob(0,:,:,:).lt.siglvl,True) ;use for noTP
prob = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
printVarSummary(var)
printVarSummary(u)
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

;----------Set res line-----------
    res_line                  = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False       ; Fill contours
    res_line@cnLinesOn        = True      ; turn off contour lines
    res_line@cnLineLabelsOn   = False      ; turn off contour line labels
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 2.5
    res_line@cnMonoLineColor  = True ;False
    res_line@cnLineColor      = "black"
    res_line@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_line@cnLevels         = (/20/)
    
;    cmap = read_colormap_file("GMT_gray") ;has 10 colors 
;    res_line@cnLinePalette  = "GMT_gray" ;has 64 colors blue green yellow,red
;    res_line@cnLinePalette  = "GMT_cool"  ;has 64 colors blue green yellow,red
;    res_line@cnLinePalette = "WhBlGrYeRe"  ;has 100 colors,start with white
;    res_line@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
;    res_line@cnMinLevelValF        = 8  ;4          ; set min contour level
;    res_line@cnMaxLevelValF        = 32 ;60 ;32        ; set max contour level
;    res_line@cnLevelSpacingF       = 8 ;        ; set contour interval
;    res_line@cnLineColors     =  (/8,6,4,2/) ;(/"lemonchiffon","lemonchiffon1","lemonchiffon2","lemonchiffon3","lemonchiffon4"/)
;    res_line@cnSpanLinePalette = True
    
    ;res_line@gsnSpreadColors      = False;True              ; turn of the full color map
    ;res_line@gsnSpreadColorStart  = 8                ; start at color 17
    ;res_line@gsnSpreadColorEnd    = 97                ; end at color 200
    res_line@lbLabelBarOn         = False   ; the label bar
    
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;--------set the attribute of polygon--------------
    lnres                      = True
    lnres@gsLineColor          = "black"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 2.0
    ;lnres@gsLineLabelString    = "TP"

    lnres2                      = True
    lnres2@gsLineColor          = "black"
    lnres2@gsLineDashPattern    = 1
    lnres2@gsLineThicknessF     = 3.5

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
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;--------set the attribute of polygon--------------
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

if(nv.eq.0) then ;preci
    minlvl = (/ 1 ,-4.0,-4.0/);-8,-8,-8,-8
    maxlvl = (/16 ,4.0 ,4.0 /);8 ,8 ,8 ,8 
    spacig = (/ 1 ,0.5 ,0.5 /);1 ,1 ,1 ,1 
end if

if(nv.eq.1.or.nv.eq.4) then ;EKE
    minlvl = (/ 2 ,-16,-16/);-8,-8,-8,-8
    maxlvl = (/50 , 16, 16/);8 ,8 ,8 ,8 
    spacig = (/ 3 , 2 , 2 /);1 ,1 ,1 ,1 
end if

if(nv.eq.2) then ;q1
    minlvl = (/0   ,-80,-80 /)
    maxlvl = (/750 , 80,  80/)
    spacig = (/ 50 , 10,  10/)
end if

if(nv.eq.3) then ;U
    minlvl = (/ 0,-4.0,-4.0/);-8,-8,-8,-8
    maxlvl = (/45, 4.0, 4.0/);8 ,8 ,8 ,8 
    spacig = (/ 3, 0.5, 0.5/);1 ,1 ,1 ,1 
end if

do nc = 0, ncase-1,1
if(nc.eq.0) then
   res@cnFillPalette = "precip3_16lev"
   ;blue and yellow and red, has 17 colors,2-18,2 is white
else
    res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
end if 
    do nm = 0,nmonth-1,1
        nn = nm +1
        np = nmonth*nc + nm
        np2= nmonth*(nc-1) + nm
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc)
        res@gsnRightString  = months(imonth(nm)) + " " +  draw_var(nv) + " (" + unit(nv) + ")" 
        ;res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")" 
        ;res@gsnRightString  = "Q1 (" + unit(nv) + ")" 
        res@gsnSpreadColors       = False              ; turn of the full color map
        if(nv.gt.10) then
        res@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
        else
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        end if
        plot(np)    = gsn_csm_contour_map_ce(wks,var(nc,nm,:,:),res)
        plot2(np)   = gsn_csm_contour(wks,u(nc,nm,:,:),res_line)
        overlay(plot(np),plot2(np))
        ;if(nc.gt.0)
        ;    plot_sig(np2)    = gsn_csm_contour(wks,prob(nc-1,nm,:,:),res_sig)
        ;    overlay(plot(np),plot_sig(np2))
        ;end if
        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        axis(np) = gsn_add_polyline( wks, plot(np), axlon, axlat(nc,nm,:), lnres2 )
        print(max(var(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nm,{lats:latn},{lonl:lonr})))
        print(max(u(nc,nm,{lats:latn},{lonl:lonr})))
        print(min(u(nc,nm,{lats:latn},{lonl:lonr})))
    end do 
end do ;for the month 

    resp@txString      = fig_title
    gsn_panel(wks,plot,(/ncase,nmonth/),resp)
end

