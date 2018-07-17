;*********************************************************************
;filename:      draw_1p5x2_clim_preci_ano.ncl
;1 calculate the climatology of the var in a cycle of two years
;2 calc the difference of pSSTclim and pSSTpCP to get the preci anomaly
;3 draw the plot
;                                                     by Ql Ren
;                                                    2018/07/04
;*********************************************************************************
begin
;area: east Asia
lats = -40 ;-20
latn = 50  ;60
lonl = 60  ;20
lonr = 280 ;220
lat_spac = 20
lon_spac = 40
lev  = (/850,500,200/)  ;can used as right string
nl = 0
siglvl = 0.1
years  = ispan(1979,2004,1)

;the path of the file
path    = "/HOME/sysu_hjkx_ys/WORKSPACE/renql/F/data"
case    = (/"pSSTclim","pSSTpEPIO","pSSTpCPIO"/)
file_in  = path + "/" + case + "/" + case + ".cam.h0"
fig_name   = path + "/fig/1p5x2_IO_sst_ano"
fig_out    ="pdf";"x11"
unit       = (/"K","K"/)   ;this can be used as right string 
var_name = (/"TS","T"/) ;this is the var_name in data file
season   = (/"(a)JJA(0)","(b)SON(0)","(c)DJF(1)","(d)MAM(1)","(e)JJA(1)"/)
nv = 0

f    = addfile(file_in(0) + ".YEAR." + years(0) + ".monthly." + var_name(0) + ".nc","r")
vars = f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
ncase  = dimsizes(case)
nsample = 13*2
nmonth  = 2*12
nseason = dimsizes(season)
nvar    = dimsizes(var_name)
nlat    = dimsizes(vars&lat) 
nlon    = dimsizes(vars&lon)
var   = new((/ncase,nsample,nmonth,nlat,nlon/),float)
var1  = new((/ncase,nsample,nseason,nlat,nlon/),float)
var_diff  = new((/ncase-1,nseason,nlat,nlon/),float)
prob      = new((/ncase-1,nseason,nlat,nlon/),float)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
    files := systemfunc("ls " + file_in(nc) + "*" + var_name(nv) + ".nc")
    do nt = 0,nsample-1,1
        print(files(2*nt:(2*nt+1)))
        f     := addfiles(files(2*nt:(2*nt+1)),"r")
        ListSetType(f,"cat")
        var(nc,nt,:,:,:) = f[:]->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
        ;var(nc,nt,:,:,:) = f[:]->$var_name(nv)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
    end do
end do

do ns = 0,nseason-1,1
    nms = ns*3 + 5
    nme = ns*3 + 7 
    var1(:,:,ns,:,:) = dim_avg_n(var(:,:,nms:nme,:,:),2) ;calc the season average
end do

var_ave  = dim_avg_n(var1,1)        ;get(ncase,nseason,nlat,nlon)
var_vari = dim_variance_n(var1,1)   ;get(ncase,nseason,nlat,nlon)
var_vari = where(var_vari.le.0.000001.and.var_vari.ge.-0.000001,0.1,var_vari)
do nc = 0,ncase-2,1
var_diff(nc,:,:,:) = var_ave(nc+1,:,:,:) - var_ave(0,:,:,:) ;get the EI anomaly
prob(nc,:,:,:)     = ttest(var_ave(nc+1,:,:,:),var_vari(nc+1,:,:,:),nsample,var_ave(0,:,:,:),var_vari(0,:,:,:),nsample,True,False)
end do

prob   = mask(prob,prob.lt.siglvl,True)
copy_VarMeta(vars,prob(0,0,:,:))
copy_VarMeta(vars,var_diff(0,0,:,:))
printVarSummary(var_diff)
printVarSummary(prob)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;the number of the panel plot
    plot     = new((/nseason*(ncase-1)/),graphic)
    plot2    = new((/nseason*(ncase-1)/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;draw the figure
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
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
    res@cnMinLevelValF        = -1.6         ; set min contour level
    res@cnMaxLevelValF        = 1.6        ; set max contour level
    res@cnLevelSpacingF       = 0.2       ; set contour interval

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = False ;True   ; the label bar
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
    res_sig@cnFillScaleF     = 0.8
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

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    res@gsnRightString  = var_name(nv) + " (" + unit(nv) + ") siglvl: " + siglvl

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.5
    resp@pmLabelBarHeightF = 0.05
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do ns = 0,nseason-1,1
do nc = 0,ncase-2,1
    np = (ncase-1)*ns+nc
    res@gsnLeftString   = season(ns) + " " + case(nc+1)
    plot(np)     = gsn_csm_contour_map_ce(wks,var_diff(nc,ns,:,:),res)
    plot2(np)    = gsn_csm_contour(wks,prob(nc,ns,:,:),res_sig)
    overlay(plot(np),plot2(np))
    print(max(var_diff(nc,ns,{lats:latn},{lonl:lonr})))
    print(min(var_diff(nc,ns,{lats:latn},{lonl:lonr})))
end do 
end do 

;    resp@txString  = case(1) + " - " + case(0)
    gsn_panel(wks,plot,(/nseason,ncase-1/),resp)
end

