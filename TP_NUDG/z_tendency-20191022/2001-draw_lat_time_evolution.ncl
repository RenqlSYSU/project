;*********************************************************************
;Draw monthly preci in every month(3*4) of one case, just simply view 
;the distribution of a variable.
;                                                     by Ql Ren
;                                                    2020/01/03
;*********************************************************************************
begin
var_option   = 2
month_option = False
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
nl = 4 ;(/925,850,700,500,200/)
fnc  = 1    ;chose which case to draw
diff = True;False ;whether draw the difference
na = 1 ;"TP","EA","NH"
month_s = 5
month_e = 8
lev = (/925,850,700,500,200/)

if(var_option.eq.0) then 
    draw_var = (/"Qd"   ,"Qeddy","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qd"   ,"Qeddy"/)  ;use to read variables
    unit     = (/"K/Day","K/Day","K/day"/)  ;this can be used as right string 
end if
if(var_option.eq.1) then 
    draw_var = (/"fQd" ,"fQeddy","A","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"f_Qd","f_Qeddy","A","total"/)  ;use to read variables
    unit     = "10~S~-10~N~" + (/"s~S~-2~N~","s~S~-2~N~","s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.2) then 
    draw_var = (/"fQd" ,"fQeddy","A","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt","dTdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
    nnv = 1;0 ;used to read dzdt induced by different forcing
    select = (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.3) then 
    draw_var = (/"sigm1"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"sigm1"/)  ;use to read variables
    unit     = (/"10~S~-5~N~J/(Pa^2*Kg)"/)  ;this can be used as right string end if
end if

;the path of the file
path    = "/home/ys17-19/renql/project/TP_NUDG"
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h"/);,"NUDG24h"/)
    case2   = (/"CTRL","NG6h-CTRL"/);, "NG24h-CTRL"/)
    nsample = 27
end if

if(case_option.eq.1) then
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    nsample = 10
end if

if(case_option.eq.2) then
    case1   = (/"F2000","NUDG"/)
    case2   = (/"F2000","NUDG-F2000"/)
    nsample = 30
end if

if(case_option.eq.2) then
    case1   = (/"F2000","NG58"/)
    case2   = (/"F2000","NG58-F2000"/)
    nsample = 30
end if

if((var_option.lt.2).or.(var_option.eq.3))
    if(month_option) then
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_month_Qd_Qeddy_A12.nc"
    jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
    u_name = "month_ave"
    else
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_daily_Qd_Qeddy_A.nc"
    jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_day_evolution_U.nc"
    u_name = "var_day"
    end if
else
    if(month_option) then
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_month_lon12.nc"
    jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
    u_name = "month_ave"
    else
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_daily.nc"
    jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_day_evolution_U.nc"
    u_name = "var_day"
    end if
end if
fig_name = path + "/z_tendency-20191022/fig/" + case1(0) + "_clim_lat_time_" + var_name(nnv) + "-" + lev(nl)  
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p"/)
fig_out    ="pdf";"x11"

;area: east Asia
area = (/"TP","EA","NH"/)
lats = (/15  ,10   , 10  /);20 ;-20
latn = (/55  ,55  , 89 /) ;55 ;60
lonl = (/70  ,110  , 10 /) ;30 ;20
lonr = (/105 ,130 ,350 /) ;150;220

f = addfile(jetfile(0),"r")
vars  = f->$u_name$(0,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
nlat   = dimsizes(vars&lat)
nlon   = dimsizes(vars&lon)
nlev   = dimsizes(lev)
nvar   = dimsizes(draw_var)
ncase  = dimsizes(case1)

if(.not.(month_option)) then 
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
var   = new((/ncase,nvar,ntime ,nlat,nlon/),double)
u     = new((/ncase,ntime,nlat,nlon/),double)
else
ntime = month_e - month_s + 1 
var   = new((/ncase,nvar,ntime,nlat,nlon/),double)
u     = new((/ncase,ntime,nlat,nlon/),double)
end if

;the number of the panel plot
    plot      = new((/nvar/),graphic)
    line1     = new((/nvar/),graphic)
    line2     = new((/nvar/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
if(month_option) then
do nc = 0, ncase-1, 1
    print(filename(nc))
    f = addfile(filename(nc),"r")
    if((var_option.lt.2).or.(var_option.eq.3))
        do nv = 0, nvar-2,1
        var(nc,nv,:,:,:)  = f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end do
    else
        var(nc,0:2,:,:,:) = f->$var_name(nnv)$(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    f = addfile(jetfile(nc),"r")
    u(nc,:,:,:) = f->month_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do 
else
do nc = 0, ncase-1, 1
    f = addfile(filename(nc),"r")
    if((var_option.lt.2).or.(var_option.eq.3))
        do nv = 0, nvar-2,1
        var(nc,nv,:,:,:)  = f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end do
    else
        var(nc,0:2,:,:,:)  = f->$var_name(nnv)$(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    f = addfile(jetfile(nc),"r")
    u(nc,:,:,:) = f->var_day(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do 
end if

var(:,nvar-1,:,:,:) = dim_sum_n(var(:,0:(nvar-2),:,:,:),1)
var_temp = dim_avg_n_Wrap(var,4)
u_temp   = dim_avg_n_Wrap(u  ,3)
printVarSummary(var_temp)
printVarSummary(u_temp)
var_temp!0 = "case"
var_temp!1 = "var"
var_temp!2 = "time"
var_temp!3 = "lat"
var_temp&case = ispan(1,ncase,1)
var_temp&var  = ispan(1,nvar,1)
var_temp&time = ispan(1,ntime,1)
u_temp!0 = "case"
u_temp!1 = "time"
u_temp&case = ispan(1,ncase,1)
u_temp&time = ispan(1,ntime,1)
printVarSummary(var_temp)
var_m = var_temp({case|:},{var|:},{lat|:},{time|:})
u_m = u_temp({case|:},{lat|:},{time|:})
printVarSummary(var_m)
printVarSummary(u_m)

lat1 = var&lat
axlat = new((/ncase,ntime/),double)
do nc = 0,ncase-1,1
    do nt = 0,ntime-1,1
        lat_ind = maxind(u_m(nc,:,nt))
        axlat(nc,nt) = lat1(lat_ind)
    end do
end do

if(var_option.eq.0) then 
    var_m = var_m*24*3600 ;convert units of Q from K/s to K/Day 
end if

if(var_option.eq.1) then 
    var_m = var_m*10^10
end if

if(var_option.eq.2) then 
if(nnv.eq.0) then
    var_m = var_m*10^3
end if
if(nnv.eq.1) then
    var_m = var_m*10^6
end if
if(nnv.eq.2) then
    var_m = var_m*10^6
end if
if(nnv.eq.3) then
    var_m = var_m*10^12
end if
end if

if(var_option.eq.3) then 
    var_m = var_m*10^5
end if
printVarSummary(var_m)

if(diff) then 
do nc = 1, ncase-1, 1
    var_m(nc,:,:,:) = var_m(nc,:,:,:) - var_m(0,:,:,:)
end do 
end if

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
    res@vpHeightF        = 0.25
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tiXAxisOn       = False
    res@tiYAxisOn       = False
    res@tiXAxisString   = "Year"
    res@tiYAxisString   = "Latitude(N)"

    res@tmXBMode        = "Explicit"
    if(month_option) then 
    res@tmXBValues      = ispan(1,(ntime),2)
    res@tmXBLabels      = (/"Jan","Mar","May","Jul","Sep","Nov"/)
    else
    res@tmXBValues      = ispan(1,(ntime),30)
    res@tmXBLabels      = (/"1 Mar","30 Mar","29 Apr","29 May","28 Jun","28 Jul","27 Aug"/)
    end if

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
    ;res@cnLineThicknessF = 0
    ;res@gsnContourZeroLineThicknessF = 4.0
    ;res@cnLineColor      = "blue"
   
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
    
    res@lbLabelBarOn         = False ;True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of polygon for westerly jet axis-------------
    lnres                      = True
    lnres@gsLineColor          = "red"
    lnres@gsLineDashPattern    = 0
    lnres@gsLineThicknessF     = 8.0

    lnres2                      = True
    lnres2@gsLineColor          = "black"
    lnres2@gsLineDashPattern    = 0
    lnres2@gsLineThicknessF     = 8.0

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
   ; res@tiMainString         = fig_title + "in" + months(m)
    ;res@tiMainFontHeightF    = 0.04
   ; res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.04
    resp@lbLabelFontHeightF = 0.015
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(var_option.eq.0) then ;Q
    minlvl = (/-2.1,-2.1,-2.1/);
    maxlvl = (/ 2.1, 2.1, 2.1/);
    spacig = (/ 0.3, 0.3, 0.3/);
end if

if(var_option.eq.1) then ;forc
    minlvl = (/-14,-14,-14,-14/);-8,-8,-8,-8,
    maxlvl = (/ 14, 14, 14, 14/); 8, 8, 8, 8,
    spacig = (/ 2 , 2 , 2 , 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.2) then ;dzdt
if(nnv.eq.0) then   ;dzdt
    minlvl = (/-3.5,-3.5,-3.5,-3.5/);
    maxlvl = (/ 3.5, 3.5, 3.5, 3.5/);
    spacig = (/ 0.5, 0.5, 0.5, 0.5/);
end if
if(nnv.eq.1) then ;dudt
    minlvl = (/-35,-35,-35,-35/);-8,-8,-8,-8,
    maxlvl = (/ 35, 35, 35, 35/); 8, 8, 8, 8,
    spacig = (/ 5 , 5 , 5 , 5 /); 1, 1, 1, 1,
end if
if(nnv.eq.2) then ;dTdt
if((fnc.eq.0)) then  ;.and.diff 
    minlvl = (/-32,-32,-32,-32/);
    maxlvl = (/ 32, 32, 32, 32/);
    spacig = (/ 4 , 4 , 4 , 4 /);
else
    minlvl = (/-8,-8,-8,-8/);-16,-16,-16,-16
    maxlvl = (/ 8, 8, 8, 8/); 16, 16, 16, 16
    spacig = (/ 1, 1, 1, 1/); 2 , 2 , 2 , 2 
end if
end if
if(nnv.eq.3) then ;dbdt
if((fnc.eq.0)) then  ;.and.diff 
    minlvl = (/-24,-24,-24,-24/);
    maxlvl = (/ 24, 24, 24, 24/);
    spacig = (/  3,  3,  3,  3/);
else
    minlvl = (/-24,-24,-24,-24/);
    maxlvl = (/ 24, 24, 24, 24/);
    spacig = (/  3,  3,  3,  3/);
end if
end if
end if

if(var_option.eq.3) then ;forc
    minlvl = (/-0.4/);-4 
    maxlvl = (/ 0.4/); 4 
    spacig = (/0.05/);0.5
end if
    
    ;res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
    ;res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nv = 0,nvar-1,1
    np = nv
   ; np = nv*nmonth + nm
    res@gsnLeftString   = pre_case(np) + ") " + case2(fnc) 
    if(var_option.eq.2) then 
    res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nnv) + ")"
    else
    res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"
    end if
    res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np) = gsn_csm_lat_time(wks,var_m(fnc,nv,:,:),res)
    line1(np)   = gsn_add_polyline( wks, plot(np),ispan(1,(ntime),1) , axlat(0,:), lnres2) ;black
    line2(np)  = gsn_add_polyline( wks, plot(np),ispan(1,(ntime),1) , axlat(1,:), lnres)  ;red
    print(max(var_m(fnc,nv,:,:)))
    print(min(var_m(fnc,nv,:,:)))
end do 

    gsn_panel(wks,plot,(/nvar,1/),resp)
end

