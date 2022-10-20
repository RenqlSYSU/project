;*************************************************************************
;read the 27 years PRECC and PRECL of three cases at East Asia(20N­45N,110-130E)
;then calc the zonal average and the climatology (ncase,nlat,nday)
;                                            by Ql_Ren
;                                           2019/01/21
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/rain_daily/rain.daily.interim."; + years(nt) + ".nc"
var_name  = (/"r"/)

fig_name  = path2 + "fig/Obs_composit_preci_lat_time"
pre_case  = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
fig_out   ="pdf"       ;"x11"
unit      = "mm/day"   ;this can be used as right string 
draw_var = (/"preci"/)  ;this is the variable used to draw, can used as center string
case1 = "ERA-Interim"
case2 = (/"More","Less","Less-More"/)

years    = ispan(1979,2010,1)
years@calendar = "standard"
nyear = dimsizes(years)

;area:   TP , east Asia
lats = (/25 , 20/)  ;-20
latn = (/40 , 45/)  ;60
lonl = (/70 ,110/)  ;20
lonr = (/105,122/)  ;130 ;220
lat_spacing = 5
na = 1

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
f  = addfile(path2 + "mdata/" + case1 + "_intEKE_TP_ts.nc","r")
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.5.and.time(:,1).le.8)

ts = dim_avg_n(f->EKEts(0:(nyear-1),select_time),1)
ts = dim_standardize_n(ts,1,0)
pslc_year = ind(ts.ge.1.2)
nslc_year = ind(ts.le.-1.2)
pnyear = dimsizes(pslc_year)
nnyear = dimsizes(nslc_year)
print("the year of more EKE is " + pnyear)
print("the year of less EKE is " + nnyear)

f    = addfile(path4 + years(0) + ".nc", "r")
vars = f->$var_name(0)$(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
ntime = dimsizes(select_time)
var   = new((/nyear ,nlat,ntime/),float)
var1  = new((/3     ,nlat,ntime/),float)
plot  = new((/ 3 /),graphic)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    f     = addfile(path4 + years(nt) + ".nc","r")
if(isleapyear(years(nt))) then ;leap year,366 day
    select = select_time + 1
else
    select = select_time
end if 
    temp := f->$var_name$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    temp := dim_avg_n_Wrap(temp,2)
    var(nt,:,:) = temp({lat|:},{time|:})
end do   ;the loop of 3 cases 
var1(2,:,:) = dim_avg_n(var,0) ;get var(ncase,nlat,ntime)
var1(0,:,:) = dim_avg_n(var(pslc_year,:,:),0) ;get var(ncase,nlat,ntime)
var1(1,:,:) = dim_avg_n(var(nslc_year,:,:),0) ;get var(ncase,nlat,ntime)
var1(0,:,:) = var1(0,:,:) - var1(2,:,:)
var1(1,:,:) = var1(1,:,:) - var1(2,:,:)
var1(2,:,:) = var1(1,:,:) - var1(0,:,:)

copy_VarMeta(temp(0,:),var1(0,:,0))
printVarSummary(var1)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name+lonl(na)+"-"+lonr(na))  
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white

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
;    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = -4          ; set min contour level
    res@cnMaxLevelValF        =  4;32        ; set max contour level
    res@cnLevelSpacingF       = 0.5        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString  = draw_var + " (" + unit +") "+ lonl(na)+" - "+lonr(na)+"E"   

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do np = 0,2,1
    res@gsnLeftString   = pre_case(np) + ") " + case2(np)
    plot(np)  = gsn_csm_lat_time(wks,var1(np,:,:),res)
    print(max(var1(np,:,:)))
    print(min(var1(np,:,:)))
end do 
    gsn_panel(wks,plot,(/3,1/),resp)
end

