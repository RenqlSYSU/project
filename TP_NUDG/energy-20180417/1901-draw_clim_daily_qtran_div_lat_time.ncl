;*********************************************************************
;filename:      draw_1p2x3_qtran_int.ncl
;input_data:    CTRL-JJA_ano_850hpa_qu.nc 
;                                                     by Ql Ren
;                                                    2018/01/18
;*********************************************************************************
begin
draw_var = (/"climate","steady","transient"/) 
var_name_u = (/"qu","ste_qu","edd_qu"/)
var_name_v = (/"qv","ste_qv","edd_qv"/)
unit     = "10~S~-5~N~kg/(s*m~S~2~N~)" ;"K"   ;this can be used as right string 

;the path of the file
path     = "/home/ys17-19/renql/project/TP_NUDG"
filename = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_JJA_qtran_int_rmvannual.nc"
gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
fig_name = path + "/energy-20180417/fig/1p3x3_clim_daily_qtran_div_lat_time"
    case   = (/" CTRL"," NG6h"," NG24h"/)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i"/)
fig_out    ="pdf";"x11"
;fig_out    ="eps";"x11"

;area:   TP , east Asia
lats = (/25 , 20/)  ;-20
latn = (/40 , 45/)  ;60
lonl = (/70 ,110/)  ;20
lonr = (/105,122/)  ;130 ;220
lat_spacing = 5
na = 1

f    = addfile(filename(0),"r")
vars:= f->$var_name_u(0)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)})
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.5.and.time(:,1).le.8)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
ntime = dimsizes(select_time)
ncase = dimsizes(case)
nvar  = dimsizes(draw_var)
var_u = new((/ntime,nlat,nlon/),float)
var   = new((/ncase,nvar,nlat,ntime/),float)

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
phis  := conform(var_u,phis,(/1,2/))

;the number of the panel plot
    plot     = new((/ncase*nvar/),graphic)

;================================================================
;read data
;====================================================================
do nf = 0 ,ncase-1, 1
f = addfile(filename(nf),"r")
do nv = 0,nvar-1,1
    print("read " + case(nf) + " " + var_name_u(nv))
    print("read " + case(nf) + " " + var_name_v(nv))
    var_u = f->$var_name_u(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    var_v = f->$var_name_v(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    divg := uv2dv_cfd ( var_u, var_v, var_u&lat, var_v&lon, 2)
    if (na.eq.0) then ;select TP area
        divg := mask(divg,(phis.ge.1500),True)
    end if
    copy_VarMeta(var_u,divg)
    divg := dim_avg_n_Wrap(divg, 2) ;get divg(ntime,nlat)
    var(nf,nv,:,:) = divg({lat|:},{time|:})*100000
end do
end do
copy_VarMeta(var_u(0,:,0),var(0,0,:,0))
printVarSummary(divg)
printVarSummary(var)

;do nc = 1, (ncase-1),1
;    var(nc,:,:,:) = var(0,:,:,:) - var(nc,:,:,:)
;    case   = (/" CTRL"," CTRL-NG6h"," CTRL-NG24h"/)
;end do

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
    ;res@tmXBMinorValues = ispan(0,(ntime-1),2)

    ;res@tmYLMode        = "Manual"
    ;res@tmYLTickStartF  = lats(na)
    ;res@tmYLTickEndF    = latn(na)
    ;res@tmYLTickSpacingF= lat_spacing
    
    res@tmXMinorGrid                   = True
    res@tmXMinorGridLineDashPattern    = 2
    res@tmXMinorGridLineColor          = "grey0"
    res@tmXMinorGridThicknessF         = 2.0
    
    res@tmXMajorGrid                   = True
    res@tmXMajorGridLineDashPattern    = 2
    res@tmXMajorGridLineColor          = "grey0"
    res@tmXMajorGridThicknessF         = 2.0
    
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
    res@lbOrientation        = "Horizontal"
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
    resp@gsnPaperOrientation = "landscape"  ;"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(na.eq.1) then
    maxlvl = (/(/24 ,24 , 8/),(/24 ,24 , 8/),(/24 ,24 , 8/)/) ;used for East Asian moist transiation
    minlvl = (/(/-24,-24,-8/),(/-24,-24,-8/),(/-24,-24,-8/)/) ;used for CTRL-NUDG
    spacig = (/(/ 3 ,3  , 1/),(/ 3 ,3  , 1/),(/ 3 ,4  , 1/)/)
else
    maxlvl = (/(/24 ,24 , 4 /),(/24 ,24 , 4 /),(/24 ,24 , 4 /)/) ;used for East Asian moist transiation
    minlvl = (/(/-24,-24,-4 /),(/-24,-24,-4 /),(/-24,-24,-4 /)/) ;used for CTRL-NUDG
    spacig = (/(/ 3 ,3  ,0.5/),(/ 3 ,3  ,0.5/),(/ 3 ,4  ,0.5/)/)
end if

do nc = 0, ncase-1,1
;if(nc.eq.0) then
;    res@cnFillPalette = "precip3_16lev"
;    ;blue and yellow and red, has 17 colors,2-18,2 is white
;else
    res@cnFillPalette = "BlueDarkRed18"
;    ;blue and red, has 18 colors,2-10 is blue
;end if 
    do nv = 0 ,nvar-1,1
        nn = nv +1
        np = nvar*nc + nv
        res@gsnLeftString   = pre_case(np) + ")" + case(nc)
        res@gsnRightString  = draw_var(nv) + " " + lonl(na)+" - "+lonr(na)+"E" 
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF      := spacig(nc,nv)        ; set contour interval
        res@cnMaxLevelValF       := maxlvl(nc,nv) ;32    ; set max contour level
        res@cnMinLevelValF       := minlvl(nc,nv)       ; set min contour level
        plot(np)   = gsn_csm_lat_time(wks,var(nc,nv,:,:),res) 
        print(max(var(nc,nv,:,:)))
        print(min(var(nc,nv,:,:)))
    end do 
end do 
    gsn_panel(wks,plot,(/ncase,nvar/),resp)
end

