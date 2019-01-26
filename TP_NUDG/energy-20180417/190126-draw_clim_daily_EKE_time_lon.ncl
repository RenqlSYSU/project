;*********************************************************************
;filename:      draw_1p2x3_qtran_int.ncl
;input_data:    CTRL-JJA_ano_850hpa_qu.nc 
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
draw_var = (/"EKE"/)
var_name = (/"EKE"/)
unit     = (/"10~S~4~N~ J"/) ;"K"   ;this can be used as right string 
nv = 0
nif = 1 ;1 mean calc difference

;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_EKE_int_rmvannual.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/1p3x3_clim_daily_EKE_time_lon"
    case   = (/" CTRL"," NG6h"," NG24h"/)
pre_case   = (/"(a","(b","(c"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

lats = (/25 , 32/)  ;-20
latn = (/32 , 40/)  ;60
lonl = (/80 , 80/)  ;20
lonr = (/122,122/)  ;130 ;220
lat_spacing = 5
lon_spacing = 20
na = 0; 1
if(na.eq.0.or.na.eq.1) then
    times = 5
    timee = 8
    tintev = 20
    tlabal = (/"1 May","21 May","10 Jun","30 Jun","20 Jul","9 Aug","29 Aug"/)
    ;tlabal = (/"1Mar","31Mar","30Apr","30May","29Jun","29Jul"/)
    ;tlabal = (/"1Jan","31Jan","2Mar","1Apr","1May","31May","30Jun","30Jul","29Aug","28Sep","28Oct","27Nov","27Dec"/)
else
    times = 6
    timee = 8
    tintev = 20
    tlabal = (/"1 Jun","21 Jun","11 Jul","31 Jul","20 Aug"/)
end if
fig_title = draw_var(nv) + " (" + unit(nv) + ") " + lats(na) + " - " + latn(na) + "N"

f    = addfile(filename(0),"r")
vars:= f->$var_name(0)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)})
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.times.and.time(:,1).le.timee)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
ntime = dimsizes(select_time)
ncase = dimsizes(case)
var   = new((/ncase,ntime,nlon/),float)

;the number of the panel plot
plot  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-1, 1
    f = addfile(filename(nf),"r")
    print("read " + case(nf) + " " + var_name(nv))
    EKE  = f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    temp = dim_avg_n_Wrap(EKE,1)
    var(nf,:,:) = temp({time|:},{lon|:})/10000
end do
copy_VarMeta(EKE(0,0,:),var(0,0,:))
printVarSummary(var)

if(nif.eq.1) then
do nc = 1, (ncase-1),1
    var(nc,:,:) = var(0,:,:) - var(nc,:,:)
    case   = (/" CTRL"," CTRL-NG6h"," CTRL-NG24h"/)
end do
end if

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   vc_thick       = 2.5
   font_index     = 22

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
    res@vpHeightF        = 0.8;0.6
    res@vpWidthF         = 0.3;0.4

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiYAxisOn       = False
    res@tiXAxisString   = "Year"
    res@tiYAxisString   = "Latitude(N)"

    res@tmYLMode        = "Explicit"
    res@tmYLValues      = ispan(0,(ntime-1),tintev)
    res@tmYLLabels      = tlabal

    res@tmXBMode        = "Manual"
    res@tmXBTickStartF  = lonl(na)
    res@tmXBTickEndF    = lonr(na)
    res@tmXBTickSpacingF= lon_spacing
    
    res@tmYMinorGrid                   = True
    res@tmYMinorGridLineDashPattern    = 2
    res@tmYMinorGridLineColor          = "grey0"
    res@tmYMinorGridThicknessF         = 2.0
    
    res@tmYMajorGrid                   = True
    res@tmYMajorGridLineDashPattern    = 2
    res@tmYMajorGridLineColor          = "grey0"
    res@tmYMajorGridThicknessF         = 2.0
    
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
    res@tmXBMajorOutwardLengthF = -0.
    res@tmYLMajorOutwardLengthF = -0.
    res@tmYRMajorOutwardLengthF = -0.
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    res@tmXBMinorOutwardLengthF = -0.
    res@tmYLMinorOutwardLengthF = -0.
    res@tmYRMinorOutwardLengthF = -0.

    res@gsnCenterString      = ""
    res@gsnLeftString      = ""
    res@gsnRightString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font
    
;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
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
    res@lbOrientation        = "Vertical";"Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2
    res@pmLabelBarOrthogonalPosF = 0.2
    res@pmLabelBarParallelPosF   = 0.5

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------

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

if(nif.eq.1) then
    if(na.eq.1) then
    minlvl = (/ 2 ,-16,-16/) ;used for Northern Asian moist transiation
    maxlvl = (/ 34, 16, 16/) ;used for CTRL-NUDG
    spacig = (/ 2 , 2 , 2 /)
    else
    minlvl = (/ 2 ,-16,-16/) ;used for Northern Asian moist transiation
    maxlvl = (/ 34, 16, 16/) ;used for CTRL-NUDG
    spacig = (/ 2 , 2 , 2 /)
    end if
else
    if(na.eq.1) then
    minlvl = (/ 2 , 2 , 2 /) ;used for Northern Asian moist transiation
    maxlvl = (/ 34, 34, 34/) ;used for CTRL-NUDG
    spacig = (/ 2 , 2 , 2 /)
    else
    minlvl = (/ 2 , 2 , 2 /) ;used for Northern Asian moist transiation
    maxlvl = (/ 34, 34, 34/) ;used for CTRL-NUDG
    spacig = (/ 2 , 2 , 2 /)
    end if
end if
    
do nc = 0, ncase-1,1
if(nif.eq.1) then
    if(nc.eq.0) then
    res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    else
    res@cnFillPalette = "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
    end if 
else
    res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
end if
    res@gsnLeftString   = pre_case(nc) + ")" + case(nc)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF      := spacig(nc)        ; set contour interval
    res@cnMaxLevelValF       := maxlvl(nc) ;32    ; set max contour level
    res@cnMinLevelValF       := minlvl(nc)       ; set min contour level
    plot(nc)   = gsn_csm_hov(wks,var(nc,:,:),res) 
    print(max(var(nc,:,:)))
    print(min(var(nc,:,:)))
end do
    resp@txString  = fig_title
    gsn_panel(wks,plot,(/1,ncase/),resp)
end

