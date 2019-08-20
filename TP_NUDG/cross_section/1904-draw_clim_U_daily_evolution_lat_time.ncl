;*************************************************************************
;read the 27 years PRECC and PRECL of three cases at East Asia(20N­45N,110-130E)
;then calc the zonal average and the climatology (ncase,nlat,nday)
;                                            by Ql_Ren
;                                           2019/01/21
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
path    = "/home/ys17-19/renql"
filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
var_name = (/"U"/)
draw_var = (/"200hPa U"/)  ;this is the variable used to draw, can used as center string
unit     = "m/s"   ;this can be used as right string 
years    = ispan(1979,2005,1)   ;the calculate time

fig_name   = path + "/project/TP_NUDG/cross_section/fig/clim_200U_daily_evolution_lat_time"
    case   = (/"(a) CTRL","(b) NG6h","(c) NG24h"/)  ;this can be used as left string
fig_out    ="pdf"       ;"x11"

;area:   TP , east Asia
lats = (/25 , 20/)  ;-20
latn = (/50 , 45/)  ;60
lonl = (/70 ,110/)  ;20
lonr = (/105,122/)  ;130 ;220
lev  = (/200/)
lat_spacing = 5
na = 0

;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.5.and.time(:,1).le.8)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
ntime = dimsizes(select_time)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var   = new((/ncase,nyear,nlat,ntime/),float)
diff  = new((/ncase,nlat,ntime/),float)
axlat = new((/ncase,ntime/),double)
plot  = new((/ncase/),graphic)
plot2 = new((/ncase/),graphic)
axis  = new((/ncase/),graphic)
axis2  = new((/ncase/),graphic)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nc = 0 , ncase-1,1
do nt = 0,nyear-1,1
    nv = 0
    f = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    u = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    temp =dim_avg_n(u,2)  ;get integrated q(time,lat)
    copy_VarMeta(u(:,:,0),temp)
    var(nc,nt,:,:) = temp({lat|:},{time|:})
end do   ;the loop of 27 years
end do   ;the loop of 3 cases 
var := dim_avg_n(var,1) ;get var(ncase,nlat,ntime)
printVarSummary(temp&time)
copy_VarMeta(temp(0,:),var(0,:,0))
copy_VarMeta(temp(0,:),diff(0,:,0))
printVarSummary(var)

lat1 = var&lat
do nc = 0,ncase-1,1
    do nt = 0,ntime-1,1
        lat_ind = maxind(var(nc,:,nt))
        axlat(nc,nt) = lat1(lat_ind)
    end do
    diff(nc,:,:) = var(nc,:,:) - var(0,:,:)
    ;case   = (/"(a) CTRL","(b) NG6h-CTRL","(c) NG24h-CTRL"/)
end do

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name+lonl(na)+"-"+lonr(na))  
    ;gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
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
    res@vpHeightF        = 0.25
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiYAxisOn       = False
    res@tiXAxisString   = "Year"
    res@tiYAxisString   = "Latitude(N)"

    res@tmXBMode        = "Explicit"
    res@tmXBValues      = ispan(0,(ntime-1),20)
    res@tmXBLabels      = (/"1 May","21 May","10 Jun","30 Jun","20 Jul","9 Aug","29 Aug"/)

    ;res@tmYLMode        = "Manual"
    ;res@tmYLTickStartF  = lats(na)
    ;res@tmYLTickEndF    = latn(na)
    ;res@tmYLTickSpacingF= lat_spacing

    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXMinorGrid                   = True
    res@tmXMinorGridLineDashPattern    = 2
    res@tmXMinorGridLineColor          = "grey0"
    res@tmXMinorGridThicknessF         = 2.0
    
    res@tmXMajorGrid                   = True
    res@tmXMajorGridLineDashPattern    = 2
    res@tmXMajorGridLineColor          = "grey0"
    res@tmXMajorGridThicknessF         = 2.0
    
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

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
   
    ;res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 10          ; set min contour level
    ;res@cnMaxLevelValF        = 58;32        ; set max contour level
    ;res@cnLevelSpacingF       = 3        ; set contour interval

    ;res@gsnSpreadColors      = True              ; turn of the full color map
    ;res@gsnSpreadColorStart  = 2                ; start at color 17
    ;res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

;----------Set res line for the zonal wind differences-----------
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
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 5
    
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF      = 3        ; set contour interval
;    res_line@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
;    res_line@cnMinLevelValF        = 8  ;4          ; set min contour level
;    res_line@cnMaxLevelValF        = 32 ;60 ;32        ; set max contour level
;    res_line@cnLevelSpacingF       = 8 ;        ; set contour interval
;    res_line@cnLineColors      =  (/8,6,4,2/) ;(/"lemonchiffon","lemonchiffon1","lemonchiffon2","lemonchiffon3","lemonchiffon4"/)
;    res_line@cnSpanLinePalette = True
    
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;--------set the attribute of polygon for westerly jet axis-------------
    lnres                      = True
    lnres@gsLineColor          = "red"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 5.0

    lnres2                      = True
    lnres2@gsLineColor          = "black"
    lnres2@gsLineDashPattern    = 0
    lnres2@gsLineThicknessF     = 5.0

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString  = "";draw_var + " (" + unit +") "+ lonl(na)+" - "+lonr(na)+"E"   
    res@gsnLeftString   = draw_var + " (" + unit +") "+ lonl(na)+" - "+lonr(na)+"E"   

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = False;True
    ;resp@lbOrientation    = "Horizontal"
    ;resp@lbLabelStride    = 2
    ;resp@lbLabelFont      = font_index              ; make labels smaller
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    ;resp@pmLabelBarWidthF  = 0.8
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@txFontHeightF = title_font 
    ;resp@txFontThicknessF = title_thick

    minlvl = (/ 0, 0, 0/) ;used for East Asian moist transiation
    maxlvl = (/45,45,45/) ;used for CTRL-NUDG
    spacig = (/ 3, 3, 3/)
    
    do nc = 0, ncase-1,1
if(nc.lt.3) then
    res@cnFillPalette := "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
else
    res@cnFillPalette := "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
end if 
        ;res@gsnLeftString   = case(nc)
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF      := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF       := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF       := minlvl(nc)       ; set min contour level
        plot(nc)   = gsn_csm_lat_time(wks,var(nc,:,:),res)
        plot2(nc)  = gsn_csm_contour(wks,diff(1,:,:),res_line)
        overlay(plot(nc),plot2(nc))
        axis(nc)   = gsn_add_polyline( wks, plot(nc),ispan(0,(ntime-1),1) , axlat(nc,:), lnres )
        axis2(nc)  = gsn_add_polyline( wks, plot(nc),ispan(0,(ntime-1),1) , axlat(1,:), lnres2 )
        print(max(var(nc,:,:)))
        print(min(var(nc,:,:)))
    end do 

    gsn_panel(wks,plot,(/ncase,1/),resp)
end

