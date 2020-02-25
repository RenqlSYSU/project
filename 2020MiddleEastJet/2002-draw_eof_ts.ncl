;***************************************************************************
;Draw the first three EOF  principal components series 
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;                                            by Ql_Ren
;                                           2020/02/23
;******************************************************************************
begin
area = (/"ME","EA","NA"/);,"GB"
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 0  ,105 , 250/);, 0   ;30 ;20
lonr = (/ 70 ,175 , 320/);, 360 ;150;220
narea = dimsizes(area)
neval = 4
years = ispan(1979,2015,1)

path   = "/home/ys17-19/renql/project/2020MiddleEastJet"
filename   = path + "/data/eof_u_6kinds-" + area + ".nc"
fig_name   = path + "/fig/eof_u_ts"
fig_out    = "pdf" ;"eps";"x11"
pre_case   = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)
plot   = new((/3/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.013
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.8
    resxy@vpHeightF = 0.2
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    
    resxy@xyDashPatterns    = (/0,0,0,0,0,0/) ;
    resxy@xyLineThicknesses = (/4.0, 4.0, 4.0, 4.0, 4.0, 4.0/) ; make second line thicker
    resxy@xyLineColors      = (/"black","red","blue","grey","darkorange3","chartreuse4"/);
    resxy@tiXAxisOn = True
    resxy@tiXAxisString = "year"
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var)); -0.4;
    ;resxy@trYMaxF = ceil(max(var)) ; 0.2 ;
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font
    resxy@tmXBLabelFontThicknessF  = tick_thick 
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 30
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 360
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(1979,2015,3)
    resxy@tmXBLabels          = ispan(1979,2015,3)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    resxy@gsnYRefLine           = 0.0
    resxy@gsnYRefLineColor      = "black"
    resxy@gsnYRefLineThicknessF = 5.0
    resxy@gsnYRefLineDashPattern= 0

    ;legend
    resxy@pmLegendDisplayMode = "NoCreate";"Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = ispan(0,5,1) 
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.05
    resxy@pmLegendHeightF  = 0.3
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.05   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.5    ; move ref vector along X-axis
    
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF   = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = "" 
    resxy@gsnRightString   = ""

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

do na = 0, narea-1, 1
    f = addfile(filename(na),"r")
    var = f->eof_ts
    resxy@trYMinF = floor(min(var)); -0.4;
    resxy@trYMaxF = ceil(max(var)) ; 0.2 ;
    do ne = 0, neval-2, 1
        ne2 = ne + 1
        resxy@gsnLeftString   = pre_case(ne) + ") EOF" + ne2 
        resxy@gsnRightString  = area(na) + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"  
        if(ne.eq.2) then
        resxy@pmLegendDisplayMode = "Always"
        else
        resxy@pmLegendDisplayMode = "NoCreate"
        end if
        plot(ne)  = gsn_csm_xy(wks,years,(/var(:,ne,:)/),resxy)
    end do
    gsn_panel(wks,plot,(/3,1/),resp)
end do
end


