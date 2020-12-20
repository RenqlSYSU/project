;*************************************************************************
;1. read the 27 years daily EKE over TP region and then calc area average
;2. read 27 years daily rain over East Asia
;3. calc nyear*nday correlation 
;                                            by Ql_Ren
;                                           2019/12/01
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
interannual = False;True ;
allyear  = False ;use 365*nyear correlation data to draw last figure
rmvannual= True ;False;if true, remove annual cycle
tna = 0
na = 0 ;0=EA; 1=NH
nv = 0
time_option = 1
fig_out   ="eps"       ;"x11"
;fig_out   ="pdf"       ;"x11"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
case1    = "ERA-Interim"
path2    = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
filename = path2 + "/mdata/" + case1 + "_Clim_month_ave_intEKE_filter.nc";
draw_var = (/"preci" ,"intEKE"    ,"uppEKE"/)  ;this is the variable used to draw, can used as center string
unit     = (/"mm/day","10~S~4~N~J","10~S~4~N~J"/)  
if(time_option.eq.0) then 
months = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
imonth = ispan(4,7,1);(/4,5,6,7/)
nm_s = 5
nm_e = 8
end if
if(time_option.eq.1) then
months = (/"JF","MA","MJ","JA"/)
imonth = ispan(2,3,1);(/4,5,6,7/)
nm_s = 5
nm_e = 8
end if
if(time_option.eq.2) then
months = (/"DJF","MAM","JJA","SON"/)
imonth = ispan(2,2,1)
nm_s = 6
nm_e = 8
end if
nmonth = dimsizes(imonth)

tarea = (/"TP","SEA","NEA"/)
tlats = (/25  ,22   ,30 /)
tlatn = (/40  ,30   ,45 /)
tlonl = (/70  ,107  ,107/)
tlonr = (/105 ,130  ,130/)
lats = (/ 20,10  /)  ;-20
latn = (/ 50,90  /)  ;60
lonl = (/90 ,0   /)  ;20
lonr = (/130,360 /)  ;130 ;220
lat_spacing = (/15,20/)
lon_spacing = (/20,30/)
siglvl = 0.05

if(nv.eq.0)
rain_option = 1 ;0 ERA, 1 CPC
if(rain_option.eq.0) then 
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/rain_daily/rain.daily.interim."; + years(nt) + ".nc"
var_name= (/"r"/)
years   = ispan(1979,2010,1)
fig_name= path2 + "fig/Obs_corr_TE_ERApreci"
else
path4   = "/home/ys17-19/data-observation/CAMS/CPC_daily_precip/precip."; + years(nt) + ".nc"
var_name= (/"precip"/)
years   = ispan(1979,2016,1)
fig_name= path2 + "fig/Obs_corr_TE_CPCpreci"
end if
end if
if(nv.ge.1) then
path4   = path2 + "/mdata/" + case1 + "_Clim_month_ave_" + draw_var(nv) + "_filter.nc"
var_name= (/"EKE"/)
years   = ispan(1979,2016,1)
fig_name= path2 + "fig/Obs_corr_TE_EKE"
end if

pre_case  = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)

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
f = addfile(gtopfile,"r")
phis = short2flt(f->z(0,{tlats(tna):tlatn(tna)},{tlonl(tna):tlonr(tna)})) ;hgt for NCEP1,unit is gpm ;z for ERA-I, unit is m2/s2
phis = phis/9.8   ;convert unit from m2/s2 to gpm

f1   = addfile(filename,"r") ;-2.5-10day
time = cd_calendar(f1->time,0)
select_time = ind(time(:,1).ge.nm_s.and.time(:,1).le.nm_e)
term = f1->EKE(0:(nyear-1),select_time,{tlats(tna):tlatn(tna)},{tlonl(tna):tlonr(tna)})
term = mask(term,phis.gt.3000,True) 
ts2 = wgt_areaave_Wrap(term,1.0,1.0,0)
delete(term)

plot  = new((/12/),graphic)
plot2 = new((/12/),graphic)
gtop  = new((/12/),graphic)

phis := short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)})) ;hgt for NCEP1,unit is gpm ;z for ERA-I, unit is m2/s2
phis = phis/9.8   ;convert unit from m2/s2 to gpm
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
if(nv.eq.0) then 
    f     = addfile(path4 + years(0) + ".nc", "r")
    vars  = f->$var_name$(0,{lats(na)-5:latn(na)+5},{lonl(na)-5:lonr(na)+5})
    nlat  = dimsizes(vars&lat)
    nlon  = dimsizes(vars&lon)
    ntime = dimsizes(select_time)
    rain  = new((/nyear,ntime,nlat,nlon/),float)
    do nt = 0,nyear-1,1
    f     = addfile(path4 + years(nt) + ".nc","r")
    if(isleapyear(years(nt))) then ;leap year,366 day
        select = select_time+1
        ;select(59:364) = select_time(59:364) + 1
    else
        select = select_time
    end if 
    rain(nt,:,:,:) = f->$var_name$(select,{lats(na)-5:latn(na)+5},{lonl(na)-5:lonr(na)+5})
    end do   ;the loop of 3 cases
    minlvl = (/-1.6/) ;-2.4 
    maxlvl = (/ 1.6/) ; 2.4 
    spacig = (/ 0.2/) ; 0.3 
else
    f    = addfile(path4, "r")
    rain = f->$var_name(0)$(0:(nyear-1),select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    nlat = dimsizes(rain&lat)
    nlon = dimsizes(rain&lon)
    rain = rain/10000
    minlvl = (/-8/) ; -4 
    maxlvl = (/ 8/) ;  4 
    spacig = (/ 1/) ; 0.5
end if
printVarSummary(rain)
prob  = new((/12,nlat,nlon/),float)
var1  = new((/12,nlat,nlon/),float)

if(rmvannual) then
ts2  = ts2  - conform(ts2 ,dim_avg_n(ts2 ,0),(/1/))
rain = rain - conform(rain,dim_avg_n(rain,0),(/1,2,3/))
fig_name = fig_name + "_rmva"
end if
if(interannual) then
fig_name = fig_name + "_interan"
end if

if(allyear) then 
nday  = dimsizes(select_time)
nsample = nyear*nday
ts  := ndtooned(ts2)
var := reshape(rain,(/nsample,nlat,nlon/))
var2  = escorc_n(var,ts,0,0)
prob2 = rtest(var2,nsample,0)
prob2 = mask(prob2,prob2.lt.siglvl,True)
copy_VarMeta(rain(0,0,:,:),var2 )
copy_VarMeta(rain(0,0,:,:),prob2)
end if

    wks = gsn_open_wks(fig_out,fig_name+lonl(na)+"-"+lonr(na)+"-"+nm_s+nm_e)  
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white

;np = -1
do lag = 1,1,1
time1 = cd_calendar(f1->time(select_time),0)
do nm = 0, nmonth-1, 1
if(time_option.eq.0) then 
    nm2_s = nm + nm_s
    nm2_e = nm + nm_s
end if
if(time_option.eq.1) then 
    nm2_s = 2*nm + nm_s
    nm2_e = 2*nm + nm_s + 1
end if
if(time_option.eq.2) then
    nm2_s = 3*nm + nm_s
    nm2_e = 3*nm + nm_s + 2
    select2 := ind(time1(:,1).ge.nm2_s.and.time1(:,1).le.nm2_e)
end if
if(time_option.eq.2.and.nm.eq.0.and.nm_s.eq.1) then 
    select2 := ind(time1(:,1).eq.12.or.time1(:,1).le.2)
else
    select2 := ind(time1(:,1).ge.nm2_s.and.time1(:,1).le.nm2_e)
end if
    nday     = dimsizes(select2)
    if(interannual) then 
    nsample = nyear
    var  = dim_avg_n(rain(:,select2,:,:),1) ;(nyear,nlat,nlon)
    ts   = dim_avg_n(ts2(:,select2),1)
    else
    nsample  = nyear*nday
    var  := reshape(rain(:,select2,:,:),(/nsample,nlat,nlon/))
    ts   := ndtooned(ts2(:,select2))
    end if

    ;term := regCoef_n(ts,var,0,0)
    if(lag.lt.0) then 
    term := regCoef_n(ts(-lag:nsample-1),var(0:nsample+lag-1,:,:),0,0)
    else
    term := regCoef_n(ts(0:nsample-lag-1),var(lag:nsample-1,:,:),0,0)
    end if
    tval := onedtond(term@tval ,(/nlat,nlon/))
    df   := onedtond(term@nptxy,(/nlat,nlon/))
    b = tval 
    b = 0.5
    var1(nm,:,:) = term 
    prob(nm,:,:) = betainc( df/(df+tval^2), df/2.0, b )

    ;var1(nm,:,:) = escorc_n(var,ts,0,0)
    ;prob(nm,:,:) = rtest(var1(nm,:,:),nsample,0)
    print("sample number is " + nsample) 
end do
var1 = var1*100000
unit_reg = unit(nv)+"/10~S~5~N~J" ;"corr";draw_var + " (" + unit +") " ;nm_s + "-" + nm_e + " " + 
prob = mask(prob,prob.lt.siglvl,True)
prob = mask(prob,var1.eq.0,False)
copy_VarMeta(rain(0,0,:,:),var1(0,:,:))
copy_VarMeta(rain(0,0,:,:),prob(0,:,:))
printVarSummary(var1)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
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
    res@gsnMajorLatSpacing = lat_spacing(na)
    res@gsnMajorLonSpacing = lon_spacing(na)
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
    res@cnLevelSpacingF       := spacig        ; set contour interval
    res@cnMaxLevelValF        := maxlvl ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl       ; set min contour level

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
    res_gtop@cnLevels         = (/3000/)
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
    res@gsnRightString  =  unit_reg ;"corr";draw_var + " (" + unit +") " ;nm_s + "-" + nm_e + " " + 

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

do nm = 0, nmonth-1, 1
    np = nm
    ;np = np + 1 
    res@gsnLeftString  = pre_case(np) + ") " + months(imonth(nm)) + " " + draw_var(nv); + " " + lag
    plot(np)  = gsn_csm_contour_map_ce(wks,var1(nm,:,:),res)
    plot2(np) = gsn_csm_contour(wks,prob(nm,:,:),res_sig)
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),plot2(np))
    overlay(plot(np),gtop(np))
    print(max(var1(nm,:,:)))
    print(min(var1(nm,:,:)))
end do

if(allyear) then 
    np = 11
    plot(np)  = gsn_csm_contour_map_ce(wks,var2,res)
    plot2(np) = gsn_csm_contour(wks,prob2,res_sig)
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),plot2(np))
    overlay(plot(np),gtop(np))
end if
    resp@txString  = "EKE leads precip" + lag + " day" 
    gsn_panel(wks,plot,(/1,2/),resp)
end do
;    gsn_panel(wks,plot,(/3,2/),resp)
end

