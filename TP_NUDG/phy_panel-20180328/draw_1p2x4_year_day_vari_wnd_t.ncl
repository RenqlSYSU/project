;*********************************************************************
;filename:      draw_1p2x4_year_day_vari_wnd_t.ncl
;input_data:    CTRL-Season_day_vari*.nc 
;
;1 read JJA data from two files season_vari(2,nlev,nlat,nlon)
;2 do the significance test
;3 draw the shaded variance and concour line of diff
;4 panel (nvari,nlev) for one variable
;                                                     by Ql Ren
;                                                    2018/03/28
;*********************************************************************************
begin
draw_var  = (/"U","V","OMEGA","T"/)
right_str = (/"U","V","omega","T"/)
unit      = (/"m~S~2~N~/s~S~2~N~","m~S~2~N~/s~S~2~N~","10~S~-3~N~ Pa~S~2~N~/s~S~2~N~","~S~o~N~C~S~2~N~"/) ;can used as right string 
var_name  = (/"season_vari"/) ;this is the var_name in data file
;var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months   = (/"DJF","MAM","JJA","SON"/)
nm = 2
nl = 2

;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
;path     = "/home/ys17-19/renql/project/TP_NUDG"
filectrl  = path + "/analysis/mdata/CTRL"   + (/"-Season_day_vari_rmvmean_","-Clim_month_ave_"/);+ draw_var(nv) + ".nc"
filenudg  = path + "/analysis/mdata/NUDG6h" + (/"-Season_day_vari_rmvmean_","-Clim_month_ave_"/);+ draw_var(nv) + ".nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/phy_panel-20180328/fig/1p2x4_year_day_vari_rmvmean_wnd_t_" + months(nm)
    case   = (/" daily"," year"/)
pre_case   = (/"(a","(b"/)
;fig_out    = "pdf";"x11"
fig_out    = "eps";"x11"

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
lev  = (/200/)
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case)
nmonth = dimsizes(months)
nvar   = dimsizes(draw_var)
nlev   = dimsizes(lev)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
var  = new((/ncase,nvar,nlat,nlon/),float)
vari = new((/ncase,nvar,nlat,nlon/),float)
prob = new((/ncase,nvar,nlat,nlon/),float)

;the number of the panel plot
    plot     = new((/ncase*nvar/),graphic)
    line     = new((/ncase*nvar/),graphic)
    gtop     = new((/ncase*nvar/),graphic)
    plot_p   = new((/ncase*nvar/),graphic)
    plot_n   = new((/ncase*nvar/),graphic)

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
   mp_thick       = 5.0
   font_index     = 22

;================================================================
;read data
;====================================================================
do nv = 0 ,nvar-1, 1
do nf = 0 ,ncase-1, 1
    f = addfile(filectrl(nf)+draw_var(nv)+".nc","r")
    var(nf,nv,:,:)  = f->$var_name$(nm,nl,:,:) ;get an ctrl (ncase,nlev,nlat,nlon)
    f = addfile(filenudg(nf)+draw_var(nv)+".nc","r")
    vari(nf,nv,:,:) = f->$var_name$(nm,nl,:,:) ;get an nudg (ncase,nlev,nlat,nlon)
end do
end do

var_temp = where(var.eq.0,var@_FillValue,var)
prob(0,:,:,:)  = ftest(var(0,:,:,:),92,vari(0,:,:,:),92,0)
prob(1,:,:,:)  = ftest(var(1,:,:,:),27,vari(1,:,:,:),27,0)
vari  = (vari - var)/var_temp*100 ;get (nudg-ctrl)/ctrl

var(:,2,:,:) = var(:,2,:,:)*1000

copy_VarMeta(var(0,0,:,:),prob(0,0,:,:))
vari = mask(vari,prob.lt.siglvl,True)
vari_p = vari
vari_n = vari
vari_p = mask(vari,vari.le.0,False)
vari_n = mask(vari,vari.lt.0,True)
printVarSummary(var)
printVarSummary(vari)
printVarSummary(prob)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"WhiteBlueGreenYellowRed") 
    ;has 254 colors,can use(/2,17,33,49,65, 81,97,113,129,145, 161,177,193,209,225, 241/)
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
    res_sigp                  = True
    res_sigp@gsnDraw          = False
    res_sigp@gsnFrame         = False
    res_sigp@cnFillOn         = True       ; Fill contours
    res_sigp@cnLinesOn        = False      ; turn off contour lines
    res_sigp@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sigp@cnInfoLabelOn    = False

    res_sigp@cnMonoFillPattern= True;False
    res_sigp@cnFillPattern    = 14
    res_sigp@cnMonoFillScale  = True
    res_sigp@cnFillScaleF     = 0.6
    res_sigp@cnMonoFillColor  = True
    res_sigp@cnFillColor      = "grey0"
    ;res_sigp@cnLinesOn        = True      ; turn off contour lines
    ;res_sigp@cnLineThicknessF = mp_thick-2
    ;res_sigp@cnLineDashPattern= 0
    ;res_sigp@cnLineColor      = "green"
    res_sigp@lbLabelBarOn     = False ;True   ; turn off the label bar
    ;res_sigp@cnLevelSelectionMode  = "Manuallevels"
    ;res_sigp@cnMinLevelValF        = 0          ; set min contour level
    ;res_sigp@cnMaxLevelValF        = 500 ;32        ; set max contour level
    ;res_sigp@cnLevelSpacingF       = 20        ; set contour interval
    
    res_sigp@gsnCenterString  = ""
    res_sigp@gsnLeftString    = " "
    res_sigp@gsnRightString   = " "
    
    res_sign    = res_sigp
    res_sign@cnFillPattern    = 11
    ;res_sigp@cnFillColor      = "green2"
    
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

minlvl = (/(/  15,  10,  0,  0/),(/  0 ,  0  ,  0  ,  0 /)/)
maxlvl = (/(/ 210, 160, 15,7.5/),(/  15,  7.5,0.75 ,0.75/)/)
spacig = (/(/  13,  10,  1,0.5/),(/   1, 0.5 ,0.05 ,0.05/)/)

do nc = 0, ncase-1,1
    do nv = 0, nvar-1,1
        nn = nv +1
        np = nvar*nc + nv
        res@gsnLeftString   = pre_case(nc) + nn + ")" + case(nc)
        res@gsnRightString  = lev + "hPa " + right_str(nv) + " (" + unit(nv) + ")"  
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc,nv)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc,nv) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc,nv)        ; set min contour level
        plot(np)    = gsn_csm_contour_map_ce(wks,var(nc,nv,:,:),res)
    
        if(all(ismissing(vari_p(nc,nv,:,:)))) then
            print("vari_p is all missing. Cannot create plot.")
        else 
            print("pertentage of (NUDG-CTRL)/CTRL*100")
            print(min(vari_p(nc,nv,{lats:latn},{lonl:lonr})))
            plot_p(np)  = gsn_csm_contour(wks,vari_p(nc,nv,:,:),res_sigp)
            overlay(plot(np),plot_p(np))
        end if

        if(all(ismissing(vari_n(nc,nv,:,:)))) then
            print("vari_n is all missing. Cannot create plot.")
        else 
            print("pertentage of (NUDG-CTRL)/CTRL*100")
            print(max(vari_n(nc,nv,{lats:latn},{lonl:lonr})))
            plot_n(np)  = gsn_csm_contour(wks,vari_n(nc,nv,:,:),res_sign)
            overlay(plot(np),plot_n(np))
        end if
        
        gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(np),gtop(np))
        line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
        print(max(var(nc,nv,{lats:latn},{lonl:lonr})))
        print(min(var(nc,nv,{lats:latn},{lonl:lonr})))
    end do 
end do
    gsn_panel(wks,plot,(/ncase,nvar/),resp)
end

