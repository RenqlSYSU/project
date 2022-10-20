;*********************************************************************
;draw temperature (color) and pressure (contour) at every time and level
;                                                     by Ql Ren
;                                                    2019/06/18
;*********************************************************************************
begin
;draw_var = (/"preci" ,"EKE"        ,"q1"            ,"U"  /)  ;this is the variable used to draw, can used as center string
;unit     = (/"mm/day","10~S~4~N~ J","10~S~-2~N~ K/S","m/s"/)    ;this can be used as right string 
lev    = (/5,10,15/) 
minlvl = (/260,215  ,350/)
maxlvl = (/305,237.5,395/)
spacig = (/  3,1.5,  3/)
pres_spac = (/50,100,50/)
nl = 1

;area: east Asia
lats = -90 ;20 ;-20
latn =  90 ;55 ;60
lonl = 0   ;30 ;20
lonr = 360 ;150;220
lat_spacing = 30
lon_spacing = 60

;the path of the file
path   = "/home/ys17-19/renql"
case1  = (/"eul_HS94_2"/)
case2  = (/"eul_HS94_2"/)
nf = 0

filename = path + "/model/" + case1(nf) + "/" + case1(nf) + ".cam.h0.*"
fig_name = path + "/project/2019cesm2_dry/test_ncl/fig/4p4x2_tc_pres_" + lev(nl) + case1(nf)
fig_title  = case2(nf) + " lev = " + lev(nl)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

;the number of the panel plot
    nplot    = 8
    plot     = new((/nplot/),graphic)
    plot2    = new((/nplot/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.015
   title_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 1.5
   font_index     = 22

;================================================================
;read data
;====================================================================
files       = systemfunc("ls " + filename) 
f           = addfiles(files,"r")
ListSetType(f,"cat")
times       = f[:]->time
ntimes = dimsizes(times)         ; number of times in the file
printVarSummary(times)

var1 = f[:]->T(:,lev(nl),:,:)
var2 = f[:]->Z3(:,lev(nl),:,:)

wks = gsn_open_wks(fig_out,fig_name)

;==============================================================
;draw figure at every time  
;==================================================================
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
;    res@gsnAddCyclic     = False
    res@vpHeightF        = 0.2
    res@vpWidthF         = 0.6

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
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

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = True;False   ; the label bar
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
;    res_line@gsnAddCyclic     = False
    res_line@cnFillOn         = False       ; Fill contours
    res_line@cnLinesOn        = True      ; turn off contour lines
    res_line@cnLineLabelsOn   = False      ; turn off contour line labels
    res_line@cnInfoLabelOn    = True
    res_line@cnLineThicknessF = 2
    res_line@cnMonoLineColor  = True ;False
    res_line@cnLineColor      = "black"
    res_line@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
    res_line@cnLevelSpacingF       = pres_spac(nl)
    
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = ""
    res_line@gsnRightString   = ""

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString  = "T(K)&Z3(m)" 

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait";"landscape" 
    ;resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.5
    resp@pmLabelBarHeightF = 0.05
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick


   res@cnFillPalette = "precip3_16lev"
   ;blue and yellow and red, has 17 colors,2-18,2 is white
   ;res@cnfillpalette = "blueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue

do np = 0,3,1
    ts = nplot*np
    te = ts + nplot - 1
    ip = 0
    do it = ts, te, 1
        res@gsnLeftString   = pre_case(ip) + ") " + times(it) 
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nl)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nl) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nl)       ; set min contour level
        
        plot(ip)    = gsn_csm_contour_map_ce(wks,var1(it,:,:),res)
        plot2(ip)   = gsn_csm_contour(wks,var2(it,:,:),res_line)
        overlay(plot(ip),plot2(ip))
        print(max(var1(it,:,:)))
        print(min(var1(it,:,:)))
        print(max(var2(it,:,:)))
        print(min(var2(it,:,:)))
        ip = ip + 1
    end do 

    resp@txString      = fig_title
    gsn_panel(wks,plot,(/4,2/),resp)
end do
end

