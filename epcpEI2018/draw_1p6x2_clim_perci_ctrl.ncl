;*********************************************************************
;filename:      draw_1p6x2_clim_preci_ctrl.ncl
;1 calculate the climatology of the var in one years
;2 draw the plot
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
lev  = (/850,500,200/)  ;can used as right string
siglvl = 0.05
years  = ispan(1979,2004,1)

;the path of the file
path    = "/home/ys17-19/renql/project"
case    = (/"GPCP"/)
file_in  = "/home/ys17-19/data-observation/GPCP-precip/precip.mon.mean.nc" 
fig_name = path + "/epcpEI2018/1p6x2_clim_preci_ctrl-" + lonl + "-" + lonr
fig_out    ="pdf";"x11"
unit       = "mm/day"   ;this can be used as right string 
draw_var = (/"preci"/)  ;this is the variable used to draw, can used as center string
;var_name = (/"PRECC","PRECL"/) ;this is the var_name in data file
var_name = (/"precip"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

f    = addfile(file_in,"r")
;f    = addfile(file_in(0) + ".YEAR." + years(0) + ".daily." + var_name(0) + ".nc","r")
;f    = addfile(file_in(0) + ".YEAR." + years(0) + ".monthly." + var_name(0) + ".nc","r")
vars = f->$var_name(0)$(0:311,{lats:latn},{lonl:lonr})
ncase  = dimsizes(case)
nsample = dimsizes(years)*12 
nmonth  = 12
nday    = 365 
nvar    = dimsizes(var_name)
nlat    = dimsizes(vars&lat) 
nlon    = dimsizes(vars&lon)
;var   = new((/ncase,nsample,nlat,nlon/),float)
;var_ave   = new((/ncase,nmonth,nlat,nlon/),float)
var_ave   = new((/nmonth,nlat,nlon/),float)

;================================================================
;read data
;====================================================================
;var = 0
;do nc = 0,ncase-1,1
;do nv = 0,nvar-1,1
;    files := systemfunc("ls " + file_in(nc) + "*" + var_name(nv) + ".nc")
;    do nt = 0,nsample-1,1
;        f     := addfile(files(nt),"r")
;       var(nc,nt,:,:,:) = var(nc,nt,:,:,:) + f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
;    end do
;end do
;end do
;var = var*3600*24*1000 ;convert unit from m/s to mm/day

;var := dim_avg_n(var,1)        ;get(ncase,nseason,nlat,nlon)

;time = cd_calendar(f->time,0)
;do nm = 0,11,1
;    select_time := ind(time(:,1).eq.nm+1)
;    var_ave(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
;end do

do nm = 0,11,1
var_ave(nm,:,:) = dim_avg_n_Wrap(var(nm:(26*12-1):12,:,:),0)   ;used for GPCP
end do
copy_VarMeta(vars,var_ave(0,:,:))
printVarSummary(var_ave)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;the number of the panel plot
    plot     = new((/nmonth/),graphic)
    plot2    = new((/nmonth/),graphic)

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
    gsn_define_colormap(wks,"precip3_16lev")
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
    res@cnMinLevelValF        = 1         ; set min contour level
    res@cnMaxLevelValF        = 16        ; set max contour level
    res@cnLevelSpacingF       = 1       ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = False   ; the label bar
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
    res@gsnRightString  = draw_var + " (" + unit + ")"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.5
    resp@pmLabelBarHeightF = 0.025
    resp@lbLabelFontHeightF = tick_font-0.01
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do ns = 0,nmonth-1,1
    res@gsnLeftString   = months(ns) + " " + case 
    plot(ns)     = gsn_csm_contour_map_ce(wks,var_ave(0,ns,:,:),res)
;    plot2(ns)    = gsn_csm_contour(wks,prob(ns,:,:),res_sig)
;    overlay(plot(ns),plot2(ns))
    print(max(var_ave(0,ns,{lats:latn},{lonl:lonr})))
    print(min(var_ave(0,ns,{lats:latn},{lonl:lonr})))
end do 

;    resp@txString  = case
    gsn_panel(wks,plot,(/6,2/),resp)
end

