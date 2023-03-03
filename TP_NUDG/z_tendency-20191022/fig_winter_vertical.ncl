;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
case_option  = 2 ;0 means NUDG experiment; 1 means noTP experiment
var_option = -4 ;0 means Qd; 1 means dzdt by 3forc; 2 advection; 3 advect; 4 high and low 
nnv = 1  ;used to read dzdt induced by different forcing

barshare = False;True ;
barvert  = False;True ;
jet_option    = True  ;20191205, control the draw of jet stream
cont_option = True ;False; used to draw dudt
diff = True ;whether draw the difference ;False;
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
fig_out    = "pdf";"pdf";"x11"

select_time_u = ispan(0,1,1)
select_time = ispan(0,1,1)
;months  = (/"JF","MJ","MJ","JA"/)
months   = (/"DJF","MAM","JJA","SON"/)
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;months  = (/"May","Jun","Jul","Aug"/)
;months  = (/"Apr","May","Jun","Jul"/)
nmonth = dimsizes(select_time)
lev   = (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,150,100/)
lev_u = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;hPa 13levels

;area: east Asia
area = (/"TP","EA","NP","NP","UP","NH","ME","EA"/)
lats = (/15  , 15 , 10 , 10 , 15 , 10 , 10 , 10 /);15 ;-210
latn = (/55  , 55 , 80 , 80 , 55 , 89 , 60 , 60 /) ;55 ;60
lonl = (/70  ,110 ,150 ,150 , 40 , 1  , 0  , 250/) ;30 ;20
lonr = (/105 ,130 ,210 ,170 , 60 ,359 , 70 , 320/) ;150;220

do na = 1,2,1 ;"TP","EA","UP","NH","CA","NP"
if(var_option.eq.-2) then 
    draw_var = (/"b1","b2","b3","vort","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"b1","b2","b3","vort","total"/)  ;this is the variable used to draw, can used as center string
    unit     = "10~S~-11~N~" + (/"s~S~-2~N~","s~S~-2~N~","s~S~-2~N~","s~S~-1~N~","s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.-4) then 
    draw_var = (/"Qte_hgh","Qte_low","Fte_hgh","Fte_low","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qeddh","Qeddl","Ah","Al","total"/)  ;this is the variable used to draw, can used as center string
    unit     = (/"K/Day","K/Day","10~S~-11~N~ s~S~-2~N~","10~S~-11~N~ s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.-3) then 
    draw_var = (/"vort_x" ,"vort_y","Tadv","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"vort_x" ,"vort_y","Tadv","total"/)  ;this is the variable used to draw, can used as center string
    unit     = "10~S~-10~N~" + (/"s~S~-2~N~","s~S~-2~N~","s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.-1) then 
    draw_var = (/"Qd" ,"Qte","Fte","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qd","Qeddy","A","total"/)  ;use to read variables
    unit     = (/"K/Day","K/Day","10~S~-11~N~ s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.0) then 
    draw_var = (/"Qd"   ,"Qte"  ,"total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qd"   ,"Qeddy","total"/)  ;use to read variables
    unit     = (/"K/Day","K/Day","K/day"/)  ;this can be used as right string 
end if

;the path of the file
path    = "/home/ys17-23/Extension2/renql/project/TP_NUDG"
if(case_option.eq.2) then
    case1   = (/"F2000","NGTR"/)
    case2   = (/"F2000","NGTR-F2000","NGTR-F2000"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 30
end if

jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
;jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_U.nc"
if(var_option.eq.-4) then 
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_season_forc_te.nc"
    contname = path+"/z_tendency-20191022/mdata/"+case1+"-Clim_dzdt_te_season_lon.nc"
    fig_name = path + "/z_tendency-20191022/fig/1p2x3_height_lat_" + var_name(0) + area(na)
end if
if(var_option.eq.-2)
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_month_advect.nc"
    fig_name = path + "/z_tendency-20191022/fig/1p2x3_height_lat_" + var_name(0) + area(na)
end if
if(var_option.le.0.and.var_option.ge.-1)
    if(nmonth.eq.1) then 
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_DJF_Qd_Qeddy_A.nc"
    else
    contname = path+"/z_tendency-20191022/mdata/"+case1+"-Clim_dzdt_season_lon.nc"
    filename = path+"/z_tendency-20191022/mdata/"+case1+"-Clim_season_Qd_Qeddy_A_lon.nc"
    ;filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_2month_Qd_Qeddy_A_lon18.nc"
    end if
    fig_name = path + "/z_tendency-20191022/fig/1p2x3_height_lat_"+var_name(0)+lonl(na)
end if

f = addfile(gtopfile,"r")
phis   := f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nlev   = dimsizes(lev)
nvar   = dimsizes(draw_var)
ncase  = dimsizes(case2)
nfile  = dimsizes(case1)
phis1  := 1013.25*(1-dim_avg_n_Wrap(phis,1)*0.0065/288.15)^5.25145
phis1  = where(phis1.le.1000, phis1, 1000)
ypoly := new((/nlat+4/),float)
xpoly := new((/nlat+4/),float)
ypoly(0:nlat-1)    = tofloat(phis1)
ypoly(nlat:nlat+3) = (/tofloat(phis1(nlat-1)),1000,1000,tofloat(phis1(0))/)
xpoly(0:nlat-1)    = tofloat(phis&lat)
xpoly(nlat:nlat+3) = (/tofloat(latn(na)),tofloat(latn(na)),tofloat(lats(na)),tofloat(lats(na))/)
delete([/phis,phis1/])

f = addfile(filename(0),"r")
if(var_option.le.0)
    if(nmonth.eq.1) then
    vars := f->$var_name(0)$({lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    else
    vars := f->$var_name(0)$(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
else
    vars := f->$var_name(nnv)$(0,0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
end if
nlev2 = dimsizes(lev_u)
var_m := new((/ncase,nvar,nmonth,nlev,nlat,nlon/),double)
var_u := new((/ncase,nmonth,nlev2,nlat/),double)
if(cont_option) then
var_c := new((/ncase,nvar,nmonth,nlev,nlat,nlon/),double)
end if

;the number of the panel plot
    plot   = new((/12/),graphic)
    line   = new((/12/),graphic)
    core   = new((/12/),graphic)
    msfp   = new((/12/),graphic)
    topo   = new((/12/),graphic)
    cont  := new((/12/),graphic)

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
do nc = 0, nfile-1, 1
    if(jet_option) then
    f = addfile(jetfile(nc),"r")
    u := f->month_ave(select_time_u,{lev_u},{lats(na):latn(na)},{lonl(na):lonr(na)})
    var_u(nc,:,:,:) = dim_avg_n_Wrap(u,3)
    end if

    f = addfile(filename(nc),"r")
    if(var_option.le.0)
    
    if(nmonth.eq.1)
        do nv = 0, nvar-2,1
        var_m(nc,nv,0,:,:,:)  = f->$var_name(nv)$({lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end do
    else
        do nv = 0, nvar-2,1
        var_m(nc,nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end do
    end if
    
    else
    
    if((var_option.eq.2).and.(nnv.eq.5)) then 
        dudt = f->$var_name(1)$(select,:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        dvdt = f->$var_name(4)$(select,:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        var_m(nc,0:(nvar-2),:,:,:,:) = uv2vr_cfd( dudt, dvdt, dudt&lat, dudt&lon, 2 )
    else
        var_m(nc,0:(nvar-2),:,:,:,:)  = f->$var_name(nnv)$(select,select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    
    end if
end do

var_m := dim_avg_n(var_m,5)
var_m(:,nvar-1,:,:,:) = dim_sum_n(var_m(:,0:(nvar-2),:,:,:),1)
copy_VarMeta(vars(:,:,0),var_m(0,0,0,:,:))
printVarSummary(var_m)

if(cont_option) then 
do nc = 0, nfile-1, 1
    fc = addfile(contname(nc),"r")
    var_c(nc,0:(nvar-2),:,:,:,:) = fc->dudt(:,select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
end do
var_c := dim_avg_n(var_c,5)*10^6
var_c(2,:,:,:,:) = var_c(1,:,:,:,:) - var_c(0,:,:,:,:)
copy_VarMeta(vars(:,:,0),var_c(0,0,0,:,:))
end if

if(jet_option) then 
    ;var_u = mask(var_u,var_u.ge.-10,True)
    copy_VarMeta(u(0,:,:,0),var_u(0,0,:,:))
    printVarSummary(var_u)
    corez = new((/nfile,nmonth/),float)
    corey = new((/nfile,nmonth/),float)
    do nc = 0, nfile-1, 1
    do nt = 0, nmonth-1, 1
        term := ndtooned(var_u(nc,nt,:,:))
        indices = ind_resolve(maxind(term),(/nlev2,nlat/))
        corez(nc,nt) = tofloat(var_u&lev(indices(0,0)))
        corey(nc,nt) = tofloat(var_u&lat(indices(0,1)))
    end do 
    end do 
    printVarSummary(indices)
    print(corez)
    print(corey)
end if

if(var_option.eq.-4) then 
    var_m(:,0:1,:,:,:) = var_m(:,0:1,:,:,:)*24*3600 ;convert units of Q from K/s to K/Day 
    var_m(:,2:3,:,:,:) = var_m(:,2:3,:,:,:)*10^11
end if

if(var_option.eq.-3) then 
    var_m = var_m*10^11
end if

if(var_option.eq.-2) then 
    var_m = var_m*10^11
end if

if(var_option.eq.-1) then 
    var_m(:,0:1,:,:,:) = var_m(:,0:1,:,:,:)*24*3600 ;convert units of Q from K/s to K/Day 
    var_m(:,2,:,:,:) = var_m(:,2,:,:,:)*10^11
end if

if(var_option.eq.0) then 
    var_m = var_m*24*3600 ;convert units of Q from K/s to K/Day 
end if

if(var_option.ge.1.and.var_option.le.4) then
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
if(nnv.eq.4) then
    var_m = var_m*10^6
end if
if(nnv.eq.5) then
    var_m = var_m*10^11
end if
end if

if(var_option.eq.5) then
if(nnv.eq.0) then  ;Qeddy
    var_m = var_m*24*3600 ;convert units of Q from K/s to K/Day 
end if
if(nnv.eq.1) then ;A
    var_m = var_m*10^11
end if
if(nnv.eq.2) then ;Sflux
    var_m = var_m*10^5
end if
end if

var_m(1,:,:,:,:) = var_m(1,:,:,:,:) - var_m(0,:,:,:,:)
var_m(2,:,:,:,:) = var_c(2,:,:,:,:)
var_u(2,:,:,:) = var_u(1,:,:,:) - var_u(0,:,:,:)

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
    res@vpHeightF        = 0.6
    res@vpWidthF         = 0.8

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
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
;    res@cnLineThicknessF = 1.0

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
    
;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    
    res@gsnPresHgtHeightLabelOn = False
    res@tiYAxisOn = False
    res@tmYRLabelsOn = False
    res@tmXBMode = "Manual"
    res@tmXBTickSpacingF   = 10
    res@tmXBTickStartF     = lats(na)
    res@tmXBTickEndF       = latn(na)

;---------set the attribute of the contour-----------------
    res_msf = True
    res_msf@gsnDraw          = False
    res_msf@gsnFrame         = False
    res_msf@cnFillOn         = False
    res_msf@cnLinesOn        = True
    res_msf@cnLineLabelsOn   = True ;False
    res_msf@cnInfoLabelOn    = False
    res_msf@cnLineThicknessF = 5.0
    res_msf@cnLineColor      = "blueviolet";"gray40""black";
    res_msf@gsnContourZeroLineThicknessF = 0
    res_msf@gsnContourNegLineDashPattern = 14
    res_msf@gsnContourPosLineDashPattern = 0
    
    ;res_msf@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    ;res_msf@cnLevels             = levels2       ; set contour interval
    res_msf@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_msf@cnLevelSpacingF      = 5.0        ; set contour interval
    res_msf@gsnPresHgtHeightLabelOn = False
    res_msf@tiYAxisOn = False
    res_msf@lbLabelBarOn     = False   ; the label bar
    res_msf@gsnCenterString  = ""
    res_msf@gsnLeftString    = " "
    res_msf@gsnRightString   = " "

;---------set the attribute of the contour-----------------
    res_line = res
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.0 ;4.0
    res_line@cnLineColor      = "black"
    res_line@gsnContourZeroLineThicknessF = 5.0
    res_line@gsnContourNegLineDashPattern = 14
    res_line@gsnContourPosLineDashPattern = 0
   
    levels = new((/10/),integer)
    levels(0) = 0
    levels(1) = 5
    levels(2:9) = ispan(10,80,10)
    res_line@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    res_line@cnLevels             = levels       ; set contour interval
    ;res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    polyres                   = True
    polyres@gsMarkerIndex     = 5        ; polymarker style
    polyres@gsMarkerSizeF     = 8.          ; polymarker size 
    polyres@gsMarkerColor     = "black"
    polyres@gsMarkerThicknessF  = 5.0 

;used to draw the black terrain
    pres                  = True
    pres@gsFillColor      = "gray"

    resp  = True    
    resp@gsnMaximize   = True
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "landscaoe";"portrait"
    resp@gsnPanelLabelBar = True
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.6
    else
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.8 ;0.05
    resp@pmLabelBarHeightF = 0.05;0.8
    end if
    else
    res@lbLabelBarOn         = True   ; the label bar
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarZone               = 0
    if(barvert) then
    res@lbOrientation        = "Vertical";
    res@pmLabelBarWidthF  = 0.1
    res@pmLabelBarHeightF = 0.6
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF  = 0.8
    res@pmLabelBarHeightF = 0.1
    res@pmLabelBarOrthogonalPosF     = 0.05 ;y direction
    res@pmLabelBarParallelPosF       = 0.5
    end if
    resp@gsnPaperOrientation = "portrait"
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(var_option.eq.-4) then
if(na.eq.1) then 
    minlvl = (/(/-2.1,-2.1,-7,-7,-7/),(/-2.1,-2.1,-7,-7,-7/),(/-7,-7,-7,-7,-7/)/)
    maxlvl = (/(/ 2.1, 2.1, 7, 7, 7/),(/ 2.1, 2.1, 7, 7, 7/),(/ 7, 7, 7, 7, 7/)/)
    spacig = (/(/ 0.3, 0.3, 1, 1, 1/),(/ 0.3, 0.3, 1, 1, 1/),(/ 1, 1, 1, 1, 1/)/)
    jet_spac = 1
    dudt_spac= 5
else
    minlvl = (/(/-2.1,-2.1,-7,-7,-7/),(/-0.7,-0.7,-2.1,-2.1,-2.1/),(/-7,-7,-7,-7,-7/)/)
    maxlvl = (/(/ 2.1, 2.1, 7, 7, 7/),(/ 0.7, 0.7, 2.1, 2.1, 2.1/),(/ 7, 7, 7, 7, 7/)/)
    spacig = (/(/ 0.3, 0.3, 1, 1, 1/),(/ 0.1, 0.1, 0.3, 0.3, 0.3/),(/ 1, 1, 1, 1, 1/)/)
    jet_spac = 2
    dudt_spac= 2
end if
end if
    
if(var_option.eq.-2) then
    minlvl = (/-14,-14,-14,-14,-14/);-8,-8,-8,-8,
    maxlvl = (/ 14, 14, 14, 14, 14/); 8, 8, 8, 8,
    spacig = (/ 2 , 2 , 2 , 2 , 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.-1) then
if(na.eq.1.or.na.eq.3) then 
    minlvl = (/(/-2.1,-2.1,-7/),(/-0.7,-0.7,-7/),(/-7,-7,-7/)/)
    maxlvl = (/(/ 2.1, 2.1, 7/),(/ 0.7, 0.7, 7/),(/ 7, 7, 7/)/)
    spacig = (/(/ 0.3, 0.3, 1/),(/ 0.1, 0.1, 1/),(/ 1, 1, 1/)/)
    jet_spac = 1
    dudt_spac= 5
else
    minlvl = (/(/-2.1,-2.1,-7/),(/-0.7,-0.7,-2.1/),(/-7,-7,-7/)/)
    maxlvl = (/(/ 2.1, 2.1, 7/),(/ 0.7, 0.7, 2.1/),(/ 7, 7, 7/)/)
    spacig = (/(/ 0.3, 0.3, 1/),(/ 0.1, 0.1, 0.3/),(/ 1, 1, 1/)/)
    jet_spac = 2
    dudt_spac= 2
end if
end if
    
if(var_option.eq.0) then 
    minlvl = (/-2.1,-2.1,-2.1/);
    maxlvl = (/ 2.1, 2.1, 2.1/);
    spacig = (/ 0.3, 0.3, 0.3/);
end if

;    res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nm = 0, nmonth-1,1
do fnc = 0,ncase-1,1
do nv = 0, nvar-2,1
    np = nv*3 + fnc ;1 var, 2 month
    ;np = nv*nmonth + nm ;1 var, 2 month
    res@gsnLeftString   = pre_case(np) + ") " + case2(fnc) + " " + months(nm) 
    if(var_option.ge.1.and.var_option.le.5) then 
    res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_var(nv) ;+ " (" + unit(nnv) + ")" ;
    else
    res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_var(nv) ;+ " (" + unit(nv) + ")"
    end if
    res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(fnc,nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(fnc,nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(fnc,nv)       ; set min contour level
    plot(np) = gsn_csm_pres_hgt(wks,var_m(fnc,nv,nm,:,:),res)
    printMinMax(var_m(fnc,nv,nm,:,:),0)
    
    if(jet_option) then 
    if(fnc.eq.2) then
    ;res_line@gsnContourZeroLineThicknessF = 8
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = jet_spac       ; set contour interval
    else
    res_line@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    res_line@cnLevels             = levels       ; set contour interval
    end if
    line(np) = gsn_csm_pres_hgt(wks,var_u(fnc,nm,:,:),res_line)
    overlay(plot(np),line(np))
    ;if(.not.(.not.(diff).and.fnc.eq.2)) then
    ;core(np) = gsn_add_polymarker(wks,plot(np),corey(fnc,nm),corez(fnc,nm),polyres)
    ;end if
    end if
    
    if(cont_option.and.fnc.lt.2) then 
    ;if((diff.and.fnc.eq.1).or.fnc.eq.2) then
    ;res_msf@cnLevelSpacingF       = dudt_spac      ; set contour interval
    ;else
    ;res_msf@cnLevelSpacingF       = 5       ; set contour interval
    ;end if
    cont(np) = gsn_csm_pres_hgt(wks,var_c(fnc,nv,nm,:,:),res_msf)
    overlay(plot(np),cont(np))
    end if
    
    topo(np) = gsn_add_polygon(wks, plot(np), xpoly, ypoly, pres)
end do 

;    resp@txString      = fig_title
    ;gsn_panel(wks,plot,(/nmonth,nvar/),resp)
end do 
    ;gsn_panel(wks,plot,(/nmonth,nvar/),resp)
    ;gsn_panel(wks,plot,(/nvar,nmonth/),resp)
    ;gsn_panel(wks,plot,(/nmonth,nvar-1/),resp)
    gsn_panel(wks,plot,(/nvar-1,3/),resp)
end do 
end do 
end

