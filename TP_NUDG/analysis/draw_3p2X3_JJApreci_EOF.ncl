;*********************************************************************************
;filename:      draw_3p2X3_JJApreci_EOF.ncl
;input_data:    JJApreci_EOFn3.nc
;Draw the first three EOF spatial patterns and principal components of Asia JJA preci data
;every panel have three figures of spatial patterns (CTRL,NUDG,DIFF) 
;and one figure of time series which has two curves (CTRL,NUDG).
;                                            by Ql_Ren
;                                           2017/12/18
;*********************************************************************************
begin
path       = "/users/yangsong3/renql/project/TP_NUDG"
filename   = path + "/analysis/mdata/" + (/"CTRL","NUDG6h"/) + "-JJApreci_EOFn3.nc"
;fig_name   = path + "/analysis/fig/3p2X3_JJApreci_daily_EOF-6h"
fig_name   = path + "/analysis/fig/3p2X3_JJApreci_EOF-6h"
fig_title  = "of Asia JJA preci(mm/day)"  ;"EOF" + neval + fig_title
      case = (/"CTRL","NUDG6h","CTRL - NUDG6h"/)
fig_out    = "pdf" ;"eps";"x11"
unit       = "mm/day"
;var_name   = (/"eof27d","eof27d_ts"/)
var_name   = (/"eof27","eof27_ts"/)

;plot area 
lats = 10
latn = 50
lonl = 110
lonr = 140

f = addfile(filename(0),"r")
ncase  = dimsizes(case)
neval  = getfilevardimsizes(f,"neval")
ntime  = getfilevardimsizes(f,"year")
;ntime  = getfilevardimsizes(f,"day")
nlat   = getfilevardimsizes(f,"lat")
nlon   = getfilevardimsizes(f,"lon")
map_var  = new((/ncase  ,neval,nlat,nlon/),float)
ts_var   = new((/ncase-1,neval,ntime/),float)
pcvar    = new((/ncase-1,neval/),float)
eval     = new((/ncase-1,neval/),float)
ts_mean  = new((/ncase-1,neval/),float)
sig      = new((/ncase-1,neval/),logical)
plot     = new((/4/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   label_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 5.0

;================================================================
;read data
;====================================================================
do nf = 0,ncase-2,1
    f = addfile(filename(nf),"r")
    map_var(nf,:,:,:) = f->$var_name(0)$       ;(ncase,neval,nlat,nlon)
    ts_var (nf,:,:)   = f->$var_name(1)$       ;(ncase,neval,nyear) 
    pcvar(nf,:)   = map_var@pcvar
    eval(nf,:)    = map_var@eval
    ts_mean(nf,:) = ts_var@ts_mean
    sig(nf,:) = eofunc_north(pcvar(nf,:),ntime,True)
end do
map_var(2,:,:,:)  = map_var(0,:,:,:) - map_var(1,:,:,:)
subtitle_color = where(sig,"red","black") 
printVarSummary(map_var)
printVarSummary(ts_var)

;timed     = cd_calendar(f->day,0)
;total_day = dimsizes(timed(:,0).eq.1979) 
;time      = fspan(1979,2006-1/total_day,ntime)

;==============================================================
;draw figure 
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "landscape"

;------------set the attribute of the map----------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 10
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;------------set the form of the coordinate-----------------------
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

;------------set the attribute of the contour----------------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -0.09       ; set min contour level
    res@cnMaxLevelValF        = 0.09       ; set max contour level
    res@cnLevelSpacingF       = 0.01        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200

    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    
;------------set the resxy for time series figure-----------
    resxy     = True
    resxy@vpXF      = 0.1
    resxy@vpYF      = 0.4
    resxy@vpWidthF  = 1.5
    resxy@vpHeightF = 1
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    
    resxy@xyDashPatterns     = (/0,0/)
    resxy@xyLineThicknesses = (/2.0, 2.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"blue","red"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = -4.0
    ;resxy@trYMaxF = 4.0
    
    resxy@tmXBLabelFontHeightF     = label_font  + 0.01
    resxy@tmXBLabelFontThicknessF  = label_thick + 0.01
    resxy@tmYRLabelFontHeightF     = label_font  + 0.01
    resxy@tmYRLabelFontThicknessF  = label_thick + 0.01

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmYLMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 3
    ;resxy@tmYLTickSpacingF   = 1.0

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    ;resxy@gsnYRefLine             = 0
    resxy@gsnYRefLineDashPatterns = (/12,12/)   
    resxy@gsnYRefLineThicknesses  = (/3.0,3.0/)
    resxy@gsnYRefLineColors       = (/"blue","red"/)
                             
    resxy@gsnCenterStringFontHeightF = subtitle_font-0.005
    resxy@gsnRightStringFontHeightF  = subtitle_font + 0.01
    resxy@gsnLeftStringFontHeightF   = subtitle_font + 0.01
    resxy@gsnCenterString  = ""
    resxy@gsnLeftString    = unit
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case(0:1)
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick + 0.01
    resxy@lgLabelFontHeightF     = label_font + 0.01
    resxy@pmLegendWidthF   = 0.2
    resxy@pmLegendHeightF  = 0.15
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.38   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.3    ; move ref vector along X-axis

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    ;res@gsnCenterStringFontHeightF= subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    resp@gsnPanelRowSpec  = True
    resp@gsnPaperOrientation = "landscape"

do ne = 0,neval-2,1
    do nc = 0,ncase-1,1
        res@gsnLeftString    = case(nc)
        if(nc.eq.2) then
            res@gsnRightString   = "" 
        else
            res@gsnRightStringFontColor   = subtitle_color(nc,ne)
            res@gsnRightString   = eval(nc,ne)+"  "+pcvar(nc,ne)+"%" 
        end if
        plot(nc) = gsn_csm_contour_map_ce(wks,map_var(nc,ne,:,:),res)
        print(max(map_var(nc,ne,{lats:latn},{lonl:lonr})))
        print(min(map_var(nc,ne,{lats:latn},{lonl:lonr})))
    end do
    
    resxy@gsnRightString   = "CTRL:" + ts_mean(0,ne) + " NUDG6h:" + ts_mean(1,ne)
    resxy@gsnYRefLine      = dim_avg_n_Wrap(ts_var(:,ne,:),1) 
    plot(3) = gsn_csm_xy(wks,ts_var&year,(/ts_var(:,ne,:)/),resxy)
    ;plot(3) = gsn_csm_xy(wks,time,(/ts_var(:,ne,:)/),resxy)

    ineval=ne+1
    resp@txString      = "EOF" + ineval + " " + fig_title 
    gsn_panel(wks,plot,(/3,1/),resp)
end do 
end

