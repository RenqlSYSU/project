;*********************************************************************
;Draw monthly preci in every month(3*4) of one case, just simply view 
;the distribution of a variable.
;                                                     by Ql Ren
;                                                    2019/07/12
;*********************************************************************************
begin
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
var_option  = 4 ;-1 forcing; 1 dudt ;4 high-low
nnv = 0 ;"dzdt","dudt","dTdt","dbdt","dvdt"

month_option = True ;False
uv_option    = False;True ;
duvdt_option = False;True ;
jet_option   = 0 ;0 = no jet, 1 = jet axis, 2 = jet contour
nl   = 1 ;(/925,850,700,500,200/)
diff = True ;whether draw the difference
na = 2 ;"TP","EA","NH"
barshare = False;True ;
barvert  = False;True

months   = (/"DJF","MAM"/);,"JJA","SON"
;months = (/"JF","MA","MJ","JA"/)
nmonth = dimsizes(months)
lev = (/925,850,700,500,200/)
area = (/"TP","EA","NH"/)
lats = (/15  ,10  , 10 /);20 ;-20
latn = (/55  ,60  , 89 /) ;55 ;60
lonl = (/70  ,30  , 0  /) ;30 ;20
lonr = (/105 ,150 ,360 /) ;150;220
month_s = 1
month_e = 8

if(var_option.eq.-6) then 
    draw_var = (/"sigm1"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"sigm1"/)  ;use to read variables
    unit     = (/"10~S~-5~N~J/(Pa^2*Kg)"/)  ;this can be used as right string end if
end if
if(var_option.eq.-4) then 
    draw_var = (/"Qeddh","Qeddl","Ah","Al","total"/)  ;this is the variable used to draw, can used as center string
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
    draw_var = (/"Qd"   ,"Qte","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"Qd"   ,"Qeddy"/)  ;use to read variables
    unit     = (/"K/Day","K/Day","K/day"/)  ;this can be used as right string 
end if
if(var_option.eq.1) then 
    draw_var = (/"Qd" ,"Qte","Fte","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt","dTdt","dbdt","dvdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)","10~S~-6~N~m/s~S~2~N~"/)  ;this can be used as right string 
    select = (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.2) then 
    draw_var = (/"b1" ,"b2","b3","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt","dTdt","dbdt","dvdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)","10~S~-6~N~m/s~S~2~N~"/)  ;this can be used as right string 
    select = (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.3) then 
    draw_var = (/"vort_x" ,"vort_y","Tadv","vort"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt","dTdt","dbdt","dvdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)","10~S~-6~N~m/s~S~2~N~"/)  ;this can be used as right string 
    select = (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.4) then 
    draw_var = (/"Qeddh","Qeddl","Ah","Al","total"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"dzdt","dudt","dTdt","dbdt","dvdt"/)  ;
    unit     = (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)","10~S~-6~N~m/s~S~2~N~"/)  ;this can be used as right string 
    select = (/0,1,2,3/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if

;the path of the file
path    = "/home/ys17-23/Extension2/renql/project/TP_NUDG"
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h"/);,"NUDG24h"/)
    case2   = (/"CTRL","NG6h-CTRL"/);, "NG24h-CTRL"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
end if
if(case_option.eq.1) then
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
end if
if(case_option.eq.2) then
    case1   = (/"F2000","NGTR"/)
    case2   = (/"F2000","NGTR-F2000"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 30
end if

if((var_option.le.0))
    if(month_option) then
    if(nmonth.eq.1) then 
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_DJF_Qd_Qeddy_A.nc"
    else
    if(var_option.ge.-1) then 
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_season_Qd_Qeddy_A_lon.nc"
    end if
    if(var_option.eq.-3) then 
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_month_advect2.nc"
    end if
    if(var_option.eq.-4) then 
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_month_forc_te.nc"
    end if
    end if
    else
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_daily_Qd_Qeddy_A.nc"
    end if
    fig_name = path + "/z_tendency-20191022/fig/1p2x3_clim_" + var_name(0) + "-" + lev(nl)  
else
    if(month_option) then
    if(nmonth.eq.1) then 
        filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_DJF.nc"
    else
        if(var_option.eq.1) then 
        filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_season_lon.nc"
        end if
        if(var_option.eq.2) then 
        filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_advect_month_season.nc"
        end if
        if(var_option.eq.3) then 
        filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_advect2_month_season.nc"
        end if
        if(var_option.eq.4) then 
        filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_dzdt_te_season_lon.nc"
        end if
    end if
    else
    filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_daily_dzdt.nc"
    end if
    fig_name = path + "/z_tendency-20191022/fig/1p2x3_clim_" + draw_var(0) + "_" + var_name(nnv) + "-" + lev(nl)  
end if
;fig_title  = "Seasonal " + draw_var(nv) + " variance in " ;+season
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p"/)
fig_out    ="pdf";"x11"

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nmonth = dimsizes(months)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nlev   = dimsizes(lev)
nvar   = dimsizes(draw_var)
ncase  = dimsizes(case1)

f = addfile(filename(0),"r")
if(.not.(month_option)) then 
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
var   = new((/ncase,nvar,ntime ,nlat,nlon/),double)
end if
var_m = new((/ncase,nvar,nmonth,nlat,nlon/),double)

if(uv_option) then 
uvar = new((/ncase,nmonth,nlat,nlon/),float)
vvar = new((/ncase,nmonth,nlat,nlon/),float)
uvari = new((/ncase,nmonth,nlat,nlon/),float)
vvari = new((/ncase,nmonth,nlat,nlon/),float)
end if
if(duvdt_option) then 
uvar = new((/ncase,nvar,nmonth,nlat,nlon/),double)
vvar = new((/ncase,nvar,nmonth,nlat,nlon/),double)
end if

;the number of the panel plot
    plot     = new((/nvar*nmonth/),graphic)
    gtop     = new((/nvar*nmonth/),graphic)
    line     = new((/nvar*nmonth/),graphic)
    vect     = new((/nvar*nmonth/),graphic)
    wind2    = new((/nvar*nmonth/),graphic)
    axis     = new((/nvar*nmonth/),graphic)
    jetl     = new((/nvar*nmonth/),graphic)

;set the area of the box
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25
lonx = (/ minlon, maxlon, maxlon, minlon, minlon /)
latx = (/ minlat, minlat, maxlat, maxlat, minlat /)

;================================================================
;read data
;====================================================================
if(jet_option.ge.1) then 
axlon = phis&lon
axlat = new((/ncase,nmonth,nlon/),double)
lat1 = phis&lat
;print(lat1)
;print(ind(lat1.eq.20))
u    = new((/ncase  ,nmonth,nlat,nlon/),float)
do nf = 0 ,ncase-1, 1
    u_f = addfile(path+"/analysis/mdata/" + case1(nf) + "-Clim_month_ave_U.nc","r")
    u(nf,:,:,:)  = u_f->season_ave(:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon):,:);
    do nm = 0, nmonth-1,1
    do inln = 0, nlon-1,1
        lat_ind = maxind(u(nf,nm,{10:48},inln))
        ;lat_ind = maxind(u(nf,nm,{20:50},inln)) + 37;ind(lat1.eq.20)
        axlat(nf,nm,inln) = lat1(lat_ind)
    end do
    end do
    printVarSummary(axlat)
end do
delete(u_f)
end if

if(uv_option) then
do nc = 0, ncase-1, 1
fu = addfile(path+"/analysis/mdata/" + case1(nc) + "-Clim_month_ave_U.nc","r")
fv = addfile(path+"/analysis/mdata/" + case1(nc) + "-Clim_month_ave_V.nc","r")
uvar(nc,:,:,:) = fu->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nseason,nlat,nlon)
vvar(nc,:,:,:) = fv->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nseason,nlat,nlon)
uvari(nc,:,:,:) = fu->season_vari(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nseason,nlat,nlon)
vvari(nc,:,:,:) = fv->season_vari(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nseason,nlat,nlon)
end do
end if

if(month_option) then
do nc = 0, ncase-1, 1
    print(filename(nc))
    f = addfile(filename(nc),"r")
    if(var_option.le.0) then
        do nv = 0, nvar-2,1
        if(nmonth.eq.1) then 
        var_m(nc,nv,0,:,:)  = f->$var_name(nv)$({lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        else
        var_m(nc,nv,:,:,:)  = f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end if
        end do
    else
        var_m(nc,0:nvar-2,:,:,:)  = f->$var_name(nnv)$(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        if(duvdt_option) then 
        uvar(nc,0:nvar-2,:,:,:)  = f->dudt(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        vvar(nc,0:nvar-2,:,:,:)  = f->dvdt(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end if 
    end if
end do 
else
do nc = 0, ncase-1, 1
    f = addfile(filename(nc),"r")
    if(var_option.le.0) then
        do nv = 0, nvar-2,1
        var(nc,nv,:,:,:)  = f->$var_name(nv)$(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end do
    else
        var(nc,0:2,:,:,:)  = f->$var_name(nnv)$(select,select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end do 
;calc ealy(56) and late(78) stage of the summer 
select_time := ind(time(:,1).ge.5.and.time(:,1).le.6)
var_m(:,:,0,:,:) = dim_avg_n_Wrap(var(:,:,select_time,:,:),2)
select_time := ind(time(:,1).ge.7.and.time(:,1).le.8)
var_m(:,:,1,:,:) = dim_avg_n_Wrap(var(:,:,select_time,:,:),2)
printVarSummary(var_m)
end if

if(var_option.eq.3) then 
var_m(:,3,:,:,:) = dim_sum_n(var_m(:,0:1,:,:,:),1)
else
var_m(:,nvar-1,:,:,:) = dim_sum_n(var_m(:,0:(nvar-2),:,:,:),1)
end if
if(duvdt_option) then 
uvar(:,nvar-1,:,:,:) = dim_sum_n(uvar(:,0:(nvar-2),:,:,:),1)
vvar(:,nvar-1,:,:,:) = dim_sum_n(vvar(:,0:(nvar-2),:,:,:),1)
uvar = uvar*100000
vvar = vvar*100000
end if

if(var_option.eq.-6) then 
    var_m = var_m*10^5
end if

if(var_option.eq.-4) then 
    var_m(:,0:1,:,:,:) = var_m(:,0:1,:,:,:)*24*3600 ;convert units of Q from K/s to K/Day 
    var_m(:,2:3,:,:,:) = var_m(:,2:3,:,:,:)*10^11
end if

if(var_option.eq.-3) then 
    var_m = var_m*10^11
end if

if(var_option.eq.-1) then 
    var_m(:,0:1,:,:,:) = var_m(:,0:1,:,:,:)*24*3600 ;convert units of Q from K/s to K/Day 
    var_m(:,2,:,:,:) = var_m(:,2,:,:,:)*10^11
end if

if(var_option.eq.0) then 
    var_m = var_m*24*3600 ;convert units of Q from K/s to K/Day 
end if

if(var_option.ge.1) then 
if(nnv.eq.0) then
    var_m = var_m*10^3
end if
if(nnv.eq.1.or.nnv.eq.4) then
    var_m = var_m*10^6
end if
if(nnv.eq.2) then
    var_m = var_m*10^6
end if
if(nnv.eq.3) then
    var_m = var_m*10^12
end if
end if
printVarSummary(var_m)

if(diff) then 
do nc = 1, ncase-1, 1
    var_m(nc,:,:,:,:) = var_m(nc,:,:,:,:) - var_m(0,:,:,:,:)
    if(uv_option) then 
    uvari = where(abs(uvari).le.0.000001,0.1,uvari)
    vvari = where(abs(vvari).le.0.000001,0.1,vvari)
    uprob = ttest(uvar(0,:,:,:),uvari(0,:,:,:),nsample,uvar(1,:,:,:),uvari(1,:,:,:),nsample,True,False)
    vprob = ttest(vvar(0,:,:,:),vvari(0,:,:,:),nsample,vvar(1,:,:,:),vvari(1,:,:,:),nsample,True,False)
    print(uprob(0,20:25,40:45))
    
    uvar(nc,:,:,:) = uvar(nc,:,:,:) - uvar(0,:,:,:)
    vvar(nc,:,:,:) = vvar(nc,:,:,:) - vvar(0,:,:,:)
    
    siglvl = 0.05
    ules = mask(uvar(1,:,:,:),((uprob.lt.siglvl).or.(vprob.lt.siglvl)),False)  ;Keep less significant values 
    vles = mask(vvar(1,:,:,:),((uprob.lt.siglvl).or.(vprob.lt.siglvl)),False)
    uvar(1,:,:,:) = mask(uvar(1,:,:,:),((uprob.lt.siglvl).or.(vprob.lt.siglvl)),True) ;Keep significant values 
    vvar(1,:,:,:) = mask(vvar(1,:,:,:),((uprob.lt.siglvl).or.(vprob.lt.siglvl)),True)
    copy_VarMeta(uvar(0,:,:,:),ules)
    copy_VarMeta(uvar(0,:,:,:),vles)
    end if
    if(duvdt_option) then 
    uvar(nc,:,:,:,:) = uvar(nc,:,:,:,:) - uvar(0,:,:,:,:)
    vvar(nc,:,:,:,:) = vvar(nc,:,:,:,:) - vvar(0,:,:,:,:)
    end if
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

;define the thickness and height of the title,tick and label
if(na.eq.2) then 
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   vcref_font     = 0.007
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 1.5
   font_index     = 22
   vc_distance  = 0.012
   vc_head_size = 0.005
else
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.015
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 3.0
   font_index     = 22
   vc_distance  = 0.02
   vc_head_size = 0.01
end if

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = 15
    res@gsnMajorLonSpacing = 30
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
    
    ;res@lbOrientation        = "Vertical";"Horizontal"
    ;res@lbLabelFontHeightF   = tick_font
    ;res@lbLabelFont          = font_index              ; make labels smaller
    ;res@lbLabelStride        = 1
    ;res@lbBoxLineThicknessF  = mp_thick-1
    ;res@pmLabelBarWidthF  = 0.03
    ;res@pmLabelBarHeightF = 0.17
    ;res@pmLabelBarZone               = 0
    ;res@pmLabelBarOrthogonalPosF     = 0.01
    ;res@pmLabelBarParallelPosF       = 0.5

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
    resvc@vcMinDistanceF             = vc_distance
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = vc_head_size
    resvc@vcLineArrowHeadMaxSizeF   = vc_head_size
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    ;resvc@vcRefMagnitudeF            = 10.0             ; define vector ref mag
    ;resvc@vcRefAnnoString2           = "10m/s"
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    resvc@vcRefAnnoOrthogonalPosF    = 0.3           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resvc2 = resvc
    resvc2@vcLineArrowColor           = "grey"

;----------Set res line-----------
    res_jet                  = True
    res_jet@gsnDraw          = False
    res_jet@gsnFrame         = False
    res_jet@gsnAddCyclic     = False
    res_jet@cnFillOn         = False       ; Fill contours
    res_jet@cnLinesOn        = True      ; turn off contour lines
    res_jet@cnLineLabelsOn   = False      ; turn off contour line labels
    res_jet@cnInfoLabelOn    = False
    res_jet@cnLineThicknessF = mp_thick + 1.5 
    res_jet@cnMonoLineColor  = True ;False
    if(uv_option.or.duvdt_option) then
    res_jet@cnLineColor      = "blueviolet";"green4";"black"
    else
    res_jet@cnLineColor      = "black"
    end if
    if(jet_option.eq.1) then 
    res_jet@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_jet@cnLevels         = (/20/)
    end if
    if(jet_option.eq.2) then 
    u = mask(u, u.lt.10 , False) 
    res_jet@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_jet@cnLevelSpacingF      = 10        ; set contour interval
    end if
    res_jet@lbLabelBarOn         = False   ; the label bar
    res_jet@gsnCenterString      = ""
    res_jet@gsnLeftString    = " "
    res_jet@gsnRightString   = " "

    res_axis                      = True
    res_axis@gsLineColor          = "black"
    res_axis@gsLineDashPattern    = 1
    res_axis@gsLineThicknessF     = mp_thick + 1.0 
    ;res_axis@gsLineLabelString    = "Westly"

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
    res@gsnStringFontHeightF = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "landscaoe";"portrait"
    resp@gsnPanelLabelBar = True
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.8
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
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.17
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.05
    res@pmLabelBarOrthogonalPosF     = 0.15 ;y direction
    res@pmLabelBarParallelPosF       = 0.5
    end if
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do fnc = 0,1,1
if(duvdt_option) then 
vcref         = (/(/5 ,2 ,2 /), \  ;the level of 925hpa    ctrl,nudg,diff
                  (/10,4 ,4 /), \  ;the level of 850hpa    ctrl,nudg,diff
                  (/7 ,2 ,2 /), \  ;the level of 700hpa    ctrl,nudg,diff
                  (/7 ,2 ,2 /), \  ;the level of 500hpa    ctrl,nudg,diff
                  (/10,4 ,4 /)/)   ;the level of 200hpa    ctrl,nudg,diff
end if

if(uv_option) then 
vcref         = (/(/10,5 ,5 /), \  ;the level of 925hpa    ctrl,nudg,diff
                  (/20,5 ,5 /), \  ;the level of 850hpa    ctrl,nudg,diff
                  (/7 ,2 ,2 /), \  ;the level of 700hpa    ctrl,nudg,diff
                  (/7 ,2 ,2 /), \  ;the level of 500hpa    ctrl,nudg,diff
                  (/40,10,10/)/)   ;the level of 200hpa    ctrl,nudg,diff
end if

if(var_option.eq.-6) then ;forc
    minlvl = (/-0.4/);-4 
    maxlvl = (/ 0.4/); 4 
    spacig = (/0.05/);0.5
end if

if(var_option.eq.-4) then
    minlvl = (/-2.1,-2.1,-14,-14,-14/);-8,-8,-8,-8,
    maxlvl = (/ 2.1, 2.1, 14, 14, 14/); 8, 8, 8, 8,
    spacig = (/ 0.3, 0.3, 2 , 2 , 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.-3) then ;forc
    minlvl = (/-14,-14,-14,-14/);-8,-8,-8,-8,
    maxlvl = (/ 14, 14, 14, 14/); 8, 8, 8, 8,
    spacig = (/ 2 , 2 , 2 , 2 /); 1, 1, 1, 1,
end if

if(var_option.eq.-1) then
    minlvl = (/-2.1,-2.1,-14,-14/);-8,-8,-8,-8,
    maxlvl = (/ 2.1, 2.1, 14, 14/); 8, 8, 8, 8,
    spacig = (/ 0.3, 0.3, 2 , 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.0) then ;Q
    minlvl = (/-2.1,-2.1,-2.1/);
    maxlvl = (/ 2.1, 2.1, 2.1/);
    spacig = (/ 0.3, 0.3, 0.3/);
end if
    
if(var_option.ge.1) then ;dzdt
if(nnv.eq.0) then   ;dzdt
    minlvl = (/-3.5,-3.5,-3.5,-3.5,-3.5/);
    maxlvl = (/ 3.5, 3.5, 3.5, 3.5, 3.5/);
    spacig = (/ 0.5, 0.5, 0.5, 0.5, 0.5/);
end if
if(nnv.eq.1) then ;dudt
    minlvl = (/-35,-35,-35,-35,-35/);-8,-8,-8,-8,
    maxlvl = (/ 35, 35, 35, 35, 35/); 8, 8, 8, 8,
    spacig = (/ 5 , 5 , 5 , 5 , 5 /); 1, 1, 1, 1,
end if
if(nnv.eq.2) then ;dTdt
if((fnc.eq.0)) then  ;.and.diff 
    minlvl = (/-28,-28,-28,-28,-28/);
    maxlvl = (/ 28, 28, 28, 28, 28/);
    spacig = (/ 4 , 4 , 4 , 4 , 4 /);
else
    minlvl = (/-7,-7,-7,-7,-7/);-16,-16,-16,-16
    maxlvl = (/ 7, 7, 7, 7, 7/); 16, 16, 16, 16
    spacig = (/ 1, 1, 1, 1, 1/); 2 , 2 , 2 , 2 
end if
end if
if(nnv.eq.3) then ;dbdt
if((fnc.eq.0)) then  ;.and.diff 
    minlvl = (/-24,-24,-24,-24,-24/);
    maxlvl = (/ 24, 24, 24, 24, 24/);
    spacig = (/  3,  3,  3,  3,  3/);
else
    minlvl = (/-24,-24,-24,-24,-24/);
    maxlvl = (/ 24, 24, 24, 24, 24/);
    spacig = (/  3,  3,  3,  3,  3/);
end if
end if
if(nnv.eq.4) then ;dvdt
    minlvl = (/-35,-35,-35,-35,-35/);-8,-8,-8,-8,
    maxlvl = (/ 35, 35, 35, 35, 35/); 8, 8, 8, 8,
    spacig = (/ 5 , 5 , 5 , 5 , 5 /); 1, 1, 1, 1,
end if
end if
    
    ;res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
    ;res@cnFillPalette = "precip3_16lev"
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nv = 0,nvar-2,1
do nm = 0,0,1
    np = nv
   ; np = nm*3 + nv
   ; np = nv*nmonth + nm
    res@gsnLeftString   = pre_case(np) + ") " + case2(fnc) + " " + months(nm) 
    if(var_option.ge.1) then 
    res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " " + var_name(nnv) + " (" + unit(nnv) + ")"
    else
    res@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"
    end if
    res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(np) = gsn_csm_contour_map_ce(wks,var_m(fnc,nv,nm,:,:),res)
    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))

    if(uv_option) then 
    resvc@vcRefMagnitudeF      = vcref(nl,fnc)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl,fnc)+"m/s"
    vect(np) = gsn_csm_vector(wks,uvar(fnc,nm,:,:),vvar(fnc,nm,:,:),resvc)
    overlay(plot(np),vect(np))
        if(diff.and.fnc.gt.0)
        resvc2@vcRefMagnitudeF      = vcref(nl,fnc)             ; define vector ref mag
        resvc2@vcRefAnnoString2     = vcref(nl,fnc)+"m/s"
        wind2(np)   = gsn_csm_vector(wks,ules(nm,:,:),vles(nm,:,:),resvc2)
        overlay(plot(np),wind2(np))
        end if
    end if
    
    if(duvdt_option) then 
    resvc@vcRefMagnitudeF      = vcref(nl,fnc)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl,fnc)+"x10~S~-5~N~m/s~S~2~N~"
    vect(np) = gsn_csm_vector(wks,uvar(fnc,nv,nm,:,:),vvar(fnc,nv,nm,:,:),resvc)
    overlay(plot(np),vect(np))
    end if
    
    if(case_option.eq.0) then 
    line(np) = gsn_add_polyline( wks, plot(np), lonx, latx, lnres )
    end if 
    print(max(var_m(fnc,nv,nm,:,:)))
    print(min(var_m(fnc,nv,nm,:,:)))
        
    if(jet_option.ge.1) then 
    if(jet_option.eq.1) then 
    axis(np) = gsn_add_polyline( wks, plot(np), axlon, axlat(nc,nm,:), res_axis )
    end if
    jetl(np) = gsn_csm_contour(wks,u(nc,nm,:,:),res_jet)
    overlay(plot(np),jetl(np))
    end if
end do 
;    gsn_panel(wks,plot,(/nmonth,1/),resp)
end do 

;    resp@txString      = fig_title
;    if(nmonth.eq.1) then
    gsn_panel(wks,plot,(/nvar-1,1/),resp)
;    else
;    gsn_panel(wks,plot,(/nmonth,3/),resp)
;    end if
end do
end

