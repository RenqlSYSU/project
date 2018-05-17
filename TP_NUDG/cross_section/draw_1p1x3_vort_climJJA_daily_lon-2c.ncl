;***************************************************************************
;filename:      draw_1p1x3_vort_daily_lon-2c.ncl
;read the 27 years U and V of three cases at East Asia(20N£­45N£¬107-130E)
;then calculate the vorticity and draw the time-lon plot to see the 
;Plateau Vortex moved eastward.
;                                            by Ql_Ren
;                                           2018/05/03
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
var_name = (/"U","V"/)
years    = ispan(1979,2005,1)   ;the calculate time

fig_name   = path + "/cross_section/fig/1p1x3_climJJAvort_daily_lon" 
    case   = (/"(a) CTRL","(b) NG6h","(c) NG6h-CTRL"/);this can be used as left string
;fig_out    ="pdf";"x11"
fig_out    ="eps";"x11"
unit       = "10~S~-6~N~ s~S~-1~N~"   ;this can be used as right string 
draw_var = (/"vort"/)  ;this is the variable used to draw, can used as center string

;area:   TP , east Asia
lats = (/25 , 20, 35,30/) ;-20
latn = (/40 , 45, 40,35/)  ;60
lonl = (/70 ,107,100,80/)  ;20
lonr = (/105,122,130,130/);130 ;220
lev  = (/700,500/)
na = 3
nl = 1

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 4.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
print(ntime)

vars  := f->$var_name(0)$(0,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
nfile = dimsizes(filein)
ncase = dimsizes(case)
var   = new((/ncase,nyear,ntime,nlon/),float)   ;define five days is a week,and the last week is seven days
plot  = new((/ncase/),graphic)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nc = 0 , nfile-1,1
do nt = 0,nyear-1,1
    nv = 0
    f = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    u = f->$var_name(nv)$(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})

    nv = 1
    f = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    v = f->$var_name(nv)$(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    
    vort  = uv2vr_cfd(u, v, v&lat, v&lon, 2) ;calculate the divergence
    var(nc,nt,:,:) = dim_avg_n_Wrap(vort,1)
end do   ;the loop of 3 cases
end do
var(2,:,:,:) = var(1,:,:,:) - var(0,:,:,:)
var := dim_avg_n_Wrap(var,1) ;get var(ncase,ntime,nlon)
var  = var*1000000  ;convert unit from m/s to mm/day
copy_VarMeta(v(:,0,:),var(0,:,:))
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name + lev(nl) + "hPa" + lats(na) + " - " + latn(na))
    gsn_define_colormap(wks,"BlueDarkRed18") 
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
    res@vpHeightF        = 0.5
    res@vpWidthF         = 0.3

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiYAxisOn       = False
    res@tiXAxisString   = "Year"
    res@tiYAxisString   = "Latitude(N)"

    ;res@tmYLMode        = "Explicit"
    ;res@tmYLValues      = ispan(0,90,10)
    ;res@tmYLLabels      = ispan(0,90,10)

; Set special resources for the time axis
    resTick                = True
    resTick@ttmFormat      = "%d %c"
    resTick@ttmAxis        = "YL"
    resTick@ttmMajorStride = 15

    ;res@tmXBMode        = "Manual"
    ;res@tmXBTickStartF  = lonl(na)
    ;res@tmXBTickEndF    = lonr(na)
    ;res@tmXBTickSpacingF= 15
    
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    ;res@tmYMajorGrid                   = True
    ;res@tmYMajorGridLineDashPattern    = 2
    ;res@tmYMajorGridLineColor          = "grey0"
    ;res@tmYMajorGridThicknessF         = 2.0
    
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
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -24          ; set min contour level
    res@cnMaxLevelValF        =  24 ;32        ; set max contour level
    res@cnLevelSpacingF       =  3        ; set contour interval

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 2

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnRightStringFontHeightF = subtitle_font
    res@gsnLeftStringFontHeightF  = subtitle_font
    res@gsnRightString  = draw_var + "(" + unit +") " 

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = True
    resp@lbOrientation    = "Vertical" ;"Horizontal"
    resp@lbLabelStride    = 2
    resp@lbLabelFont      = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = tick_font              ; make labels smaller
    resp@lbBoxLineThicknessF  = mp_thick-1
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    do nc = 0, ncase-1,1
        res@gsnLeftString   = case(nc)
        time_axis_labels( u&time,res,resTick )
        plot(nc)     = gsn_csm_hov(wks,var(nc,:,:),res)
        print(max(var(nc,:,:)))
        print(min(var(nc,:,:)))
    end do 

    gsn_panel(wks,plot,(/1,3/),resp)
end


