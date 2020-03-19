;*************************************************************************
;read the 27 years PRECC and PRECL of three cases at East Asia(20N­45N,110-130E)
;then calc the zonal average and the climatology (ncase,nlat,nday)
;                                            by Ql_Ren
;                                           2019/01/21
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
data_option = 1 ;0 ERA, 1 CPC
selecmod = 1
allyear  = False
;area:    east Asia
lats = (/ 20/)  ;-20
latn = (/ 50/)  ;60
lonl = (/90 /)  ;20
lonr = (/130/)  ;130 ;220
lat_spacing = 15
lon_spacing = 20
na = 0
nm_s = 5
nm_e = 6
siglvl = 0.05

gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
path2    = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
case1 = "ERA-Interim"
if(data_option.eq.0) then 
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/rain_daily/rain.daily.interim."; + years(nt) + ".nc"
var_name  = (/"r"/)
else
path4   = "/home/ys17-19/data-observation/CPC_daily_precip/precip."; + years(nt) + ".nc"
var_name  = (/"precip"/)
end if

fig_name  = path2 + "fig/Obs_corr_TE_preci"
pre_case  = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
fig_out   ="pdf"       ;"x11"
unit      = "mm/day"   ;this can be used as right string 
draw_var = (/"preci"/)  ;this is the variable used to draw, can used as center string

years    = ispan(1979,2016,1)
years@calendar = "standard"
nyear = dimsizes(years)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.015
   title_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 5.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f1   = addfile(path2 + "mdata/" + case1 + "_intEKE_TP_ts.nc","r") ;-2.5-10day
time = cd_calendar(f1->time,0)
select_time = ind(time(:,1).ge.nm_s.and.time(:,1).le.nm_e)
ts2 = f1->EKEts(0:(nyear-1),select_time)

f    = addfile(path4 + years(0) + ".nc", "r")
vars = f->$var_name(0)$(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
ntime = dimsizes(select_time)
rain  = new((/nyear,ntime,nlat,nlon/),float)
var1  = new((/12,nlat,nlon/),float)
prob  = new((/12,nlat,nlon/),float)
plot  = new((/12/),graphic)
plot2 = new((/12/),graphic)
gtop  = new((/12/),graphic)

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)})) ;hgt for NCEP1,unit is gpm ;z for ERA-I, unit is m2/s2
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    f     = addfile(path4 + years(nt) + ".nc","r")
if(isleapyear(years(nt))) then ;leap year,366 day
    select = select_time+1
    ;select(59:364) = select_time(59:364) + 1
else
    select = select_time
end if 
    rain(nt,:,:,:) = f->$var_name$(select,{lats(na):latn(na)},{lonl(na):lonr(na)})
end do   ;the loop of 3 cases

if(allyear) then 
nday  = dimsizes(select_time)
nsample = nyear*nday
ts  := ndtooned(ts2)
var := reshape(rain,(/nsample,nlat,nlon/))
var2  = escorc_n(var,ts,0,0)
prob2 = rtest(var2,nsample,0)
prob2 = mask(prob2,prob2.lt.siglvl,True)
copy_VarMeta(vars,var2 )
copy_VarMeta(vars,prob2)
end if

time1 = cd_calendar(f1->time(select_time),0)
do nm = 0, (nm_e-nm_s), 1
    nm2_s = nm + nm_s
    nm2_e = nm + nm_s 
    select2 := ind(time1(:,1).ge.nm2_s.and.time1(:,1).le.nm2_e)
    nday     = dimsizes(select2)
    nsample  = nyear*nday
    var  := reshape(rain(:,select2,:,:),(/nsample,nlat,nlon/))
    ts   := ndtooned(ts2(:,select2))

    term := regCoef_n(ts,var,0,0)
    tval := onedtond(term@tval ,(/nlat,nlon/))
    df   := onedtond(term@nptxy,(/nlat,nlon/))
    b = tval 
    b = 0.5
    var1(nm,:,:) = term 
    prob(nm,:,:) = betainc( df/(df+tval^2), df/2.0, b )

    ;var1(nm,:,:) = escorc_n(var,ts,0,0)
    ;prob(nm,:,:) = rtest(var1(nm,:,:),nsample,0)
    print("sample number is " + nyear*nday)
end do
var1 = var1*100000
prob = mask(prob,prob.lt.siglvl,True)
prob = mask(prob,var1.eq.0,False)
copy_VarMeta(vars,var1(0,:,:))
copy_VarMeta(vars,prob(0,:,:))
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

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
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
    
    res@tmXMinorGrid                   = False
    res@tmXMinorGridLineDashPattern    = 2
    res@tmXMinorGridLineColor          = "grey0"
    res@tmXMinorGridThicknessF         = 2.0
    
    res@tmXMajorGrid                   = False
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
    res@cnMinLevelValF        = -0.8          ; set min contour level
    res@cnMaxLevelValF        =  0.8;32        ; set max contour level
    res@cnLevelSpacingF       = 0.1        ; set contour interval

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = False ;True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

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
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.05
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
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@gsnRightString  =  "mm/day/10~S~5~N~J" ;"corr";draw_var + " (" + unit +") " ;nm_s + "-" + nm_e + " " + 

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.04
    resp@lbLabelFontHeightF = 0.02
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do nm = 0, (nm_e-nm_s), 1
    res@gsnLeftString  = pre_case(nm) + ") " + months((nm+nm_s-1))
    plot(nm)  = gsn_csm_contour_map_ce(wks,var1(nm,:,:),res)
    plot2(nm) = gsn_csm_contour(wks,prob(nm,:,:),res_sig)
    gtop(nm) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(nm),plot2(nm))
    overlay(plot(nm),gtop(nm))
    print(max(var1(nm,:,:)))
    print(min(var1(nm,:,:)))
end do

if(allyear) then 
    nm = 11
    plot(nm)  = gsn_csm_contour_map_ce(wks,var2,res)
    plot2(nm) = gsn_csm_contour(wks,prob2,res_sig)
    gtop(nm) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(nm),plot2(nm))
    overlay(plot(nm),gtop(nm))
end if
    gsn_panel(wks,plot,(/1,2/),resp)
end

