;*********************************************************************
;filename:      draw_6p5x3_stream_vort_daily-1c.ncl
;input_data:    AMIP_C5PM_TP_NUDG.cam.h1.YEAR.1998.daily.U.nc 
;Draw the contour map of monthly variance of 27 years 
;And there are three cases,1AMIP-CTRL,2TP-NUDG,3 is the difference between the above two
;the third case will do the significance testing, the dot mean the difference is significance
;The data(3variables,3cases,12month,3levs) are all save in one file
;                                                     by Ql Ren
;                                                    2018/01/08
;*********************************************************************************
begin
;the path of the file
years    = ispan(1979,2005,1)   ;the calculate time
var_name = (/"U","V"/)
lev      = (/850,500,200/)
case     = (/"CTRL","NG24h","NG6h"/)
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nt = 19
nc = 0
nl = 1
nrow    = 5  ;used for layout of panel
ncolumn = 2
npanel  = 9

path    = "/users/yangsong3/renql/project/TP_NUDG"
filename = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                   "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                   "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/6p5x3_stream_vort_daily-" + years(nt) + case(nc) + lev(nl)
fig_out  = "pdf";"x11"
unit     = "10~S~-6~N~ s~S~-1~N~"
fig_title= years(nt) + " " + lev(nl) + "hPa vort(" + unit + ") in " + case(nc)

;area: east Asia
lats = 15     
latn = 45      
lonl = 60  
lonr = 130
siglvl = 0.05

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats:latn},{lonl:lonr})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

nv = 0
f    = addfile(filename(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
time        = cd_calendar(f->time,0)
printVarSummary(time)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
u = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})

nv = 1
f = addfile(filename(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
v = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})

vort = uv2vr_cfd(u, v, v&lat, v&lon, 2)*1000000 ;calculate the vorticity
copy_VarMeta(u,vort)

;the number of the panel plot
    plot   = new((/ncolumn*nrow/),graphic)
    stream = new((/ncolumn*nrow/),graphic)
    line   = new((/ncolumn*nrow/),graphic)
    gtop   = new((/ncolumn*nrow/),graphic)

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
   title_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 4.0
   font_index     = 22

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
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

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 20
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
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False  ; True
   
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
    ;res@cnLevels             = ispan(0,60,10)
    ;res@cnFillColors         = (/-1,20,47,57,94,127,152/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -40          ; set min contour level
    res@cnMaxLevelValF        = 40 ;32        ; set max contour level
    res@cnLevelSpacingF       = 5        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 65                ; end at color 200
    
    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "vertical"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbBoxLineThicknessF  = mp_thick-1

;----------Set res st-----------
    res_st                  = True       ; use plot options
    res_st@gsnDraw          = False
    res_st@gsnFrame         = False
    res_st@gsnLeftString    = " "
    res_st@gsnRightString   = " "
    res_st@stMinArrowSpacingF = 0.008             ; arrow spacing.
    res_st@stArrowLengthF     = 0.008             ; arrow length
    res_st@stArrowStride      = 3
    res_st@stLineColor        = "black"
    res_st@stLineThicknessF   = mp_thick - 2


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
    res@gsnRightString   = " " ;draw_var + "(" + unit + ")"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.5
    resp@pmLabelBarHeightF = 0.05
    resp@lbOrientation      = "horizontal" ;"vertical"
    resp@lbLabelFontHeightF = tick_font-0.02
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFont = font_index
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

nday = (/(/0,1,2,3,4/),(/5,6,7,8,9/)/)

do np = 0,npanel-1,1
    nc = -1
    do nco = 0,ncolumn-1,1
    do nr  = 0,nrow-1,1
        nc = nc + 1
        nd = np*nrow*ncolumn + nday(nco,nr)
        res@gsnLeftString = time(select_time(nd),2) + " " + months(floattoint(time(select_time(nd),1))-1) 
        plot(nc) = gsn_csm_contour_map_ce(wks,vort(nd,:,:),res)
        stream(nc) = gsn_csm_streamline(wks,u(nd,:,:),v(nd,:,:),res_st)
        overlay(plot(nc),stream(nc))
        line(nc) = gsn_add_polyline( wks, plot(nc), lonx, latx, lnres )
        print(max(vort(nd,{lats:latn},{lonl:lonr})))
        print(min(vort(nd,{lats:latn},{lonl:lonr})))
        gtop(nc) = gsn_csm_contour(wks,phis,res_gtop)
        overlay(plot(nc),gtop(nc))
    end do 
    end do 

    resp@txString      = fig_title
    gsn_panel(wks,plot,(/nrow,ncolumn/),resp)
end do
end

