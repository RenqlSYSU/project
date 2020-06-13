;***************************************************************************
;Draw the regression pattern based on the first three eof_ts 
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;                                            by Ql_Ren
;                                           2020/02/25
;******************************************************************************
begin
lev  = (/850,500,200,100/)
nl = 0
global = True
bar_share  = False;True ;
var_option = 1  ;0 means Qd; 1 means f_Qd; 2 means dzdt; 3 means duTdx; 4 means dzdt induced by b1
do var_option = -1,-3,1
if(var_option.eq.-10) then 
    draw_var := (/"uwnd"/)  ;this is the variable used to draw, can used as center string
    var_name := (/"season_ave"/)  ;this is the variable used to draw, can used as center string
    unit     := (/"m/s"/)  ;this can be used as right string 
end if
if(var_option.eq.-3) then 
    draw_var := (/"vort","Tadvect"/)  ;,"b3y" this is the variable used to draw, can used as center string
    var_name := (/"b1x","b1y"/)  ;,"b3y" use to read variables
    unit     := (/"10~S~-11~N~s~S~-2~N~","K/day"/)  ;this can be used as right string 
end if
if(var_option.eq.-2) then 
    draw_var := (/"Qte_low","Fte_low","Qte_hgh","Fte_hgh"/)  ;this is the variable used to draw, can used as center string
    var_name := (/"Qeddl","Al","Qeddh","Ah"/)  ;use to read variables
    unit     := (/"K/Day","10~S~-11~N~s~S~-2~N~","K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.-1) then 
    draw_var := (/"Qd","Qte","Fte"/)  ;this is the variable used to draw, can used as center string
    var_name := (/"Qd"   ,"Qeddy","A"/)  ;use to read variables
    unit     := (/"K/Day","K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
end if
if(var_option.eq.0) then 
    draw_var := (/"Qd"   ,"Qeddy","total"/)  ;this is the variable used to draw, can used as center string
    var_name := (/"Qd"   ,"Qeddy"/)  ;use to read variables
    unit     := (/"K/Day","K/Day","K/day"/)  ;this can be used as right string 
end if
if(var_option.ge.1) then
var_name := (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
unit     := (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
nnv = 1 ;used to read dzdt induced by different forcing
if(var_option.eq.1) then 
    draw_var := (/"Qd","Qte","Fte"/)  ;this is the variable used to draw, can used as center string
    select := (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.2) then 
    draw_var := (/"Qte_low","Fte_low","Qte_hgh","Fte_hgh"/)  ;this is the variable used to draw, can used as center string
    select := (/1,3,0,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
if(var_option.eq.3) then 
    draw_var := (/"vort","Tadvect"/)  ;,"advect","b3y" this is the variable used to draw, can used as center string
    select := (/0,1/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
end if
end if
nvar = dimsizes(draw_var)

area = (/"ME","EA","NA"/);,"GB"
lats = (/ 10 , 10 , 10 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 90 , 90 , 90 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 0  , 0  ,0   /);, 0   ;30 ;20  0  ,105 , 250
lonr = (/ 360,360 ,360 /);, 360 ;150;220 70 ,175 , 320
lnld = (/ 0  ,105 , 250/) ;draw vertical line
lnrd = (/ 70 ,175 , 320/) ;
lonx = (/ lnld, lnrd, lnrd, lnld, lnld /)
latx = (/ lats, lats, latn, latn, lats /)
lon_spacing = 30
lat_spacing = 20
narea = dimsizes(area)
neval = 4
siglvl = 0.05
years = ispan(1979,2015,1)

case   = "ERA-Interim"
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"
filename0  = path + "data/eof_u_6kinds-" + area + ".nc"
if(var_option.eq.-10)
    filename = path + "data/" + case + "_DJF_" + draw_var(0) + "-year.nc"
    fig_name   = path + "fig/eof_regression_hor-" + draw_var(0)+lev(nl)
end if
if(var_option.eq.-3)
    filename = path + "data/" + case + "_DJF_forc_in2-year.nc"
    fig_name = path + "fig/eof_regression_hor-" + var_name(0)+lev(nl)
end if
if(var_option.eq.-2)
    filename = path + "data/" + case + "_DJF_forc_te-year.nc"
    fig_name = path + "fig/eof_regression_hor-" + var_name(0)+lev(nl)
end if
if(var_option.le.0.and.var_option.ge.-1)
    filename = path + "data/" + case + "_DJF_forc_ex-year.nc"
    fig_name = path + "fig/eof_regression_hor-" + var_name(0)+lev(nl)
end if
if(var_option.eq.1)
    filename = path + "data/" + case + "_DJF_dzdt_ex-year.nc"
    fig_name = path + "fig/eof_regression_hor-" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.2)
    filename = path + "data/" + case + "_DJF_dzdt_te-year.nc"
    fig_name = path + "fig/eof_regression_hor-" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.3)
    filename = path + "data/" + case + "_DJF_dzdt_in2-year.nc"
    fig_name = path + "fig/eof_regression_hor-" + var_name(nnv) + "_" + draw_var(0)
end if
fig_out    = "pdf" ;"eps";"x11"
pre_case   = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i", "(j","(k","(l", "(m","(n","(o", "(p","(q","(r"/)
plot   = new((/18/),graphic)
line   = new((/18/),graphic)
regl   = new((/18/),graphic)
sigl   = new((/18/),graphic)
line1  = new((/12/),graphic)

;===============================================
;draw the figure
;=============================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.02
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font

;--------set the form of the coordinate-----------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False

    res@lbOrientation        = "Vertical";"Horizontal"
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.15
    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5

;-------set the attribute of the map-------------------
    if(global) then 
    res@gsnAddCyclic     = True 
    res@mpCenterLonF     = 150.
    res@mpMaxLatF        = latn(0) 
    res@mpMinLatF        = lats(0)
    res@mpMaxLonF        = 330 ;lonr(na)
    res@mpMinLonF        = -30 ;lonl(na)
    else
    res@gsnAddCyclic     = False
    res@mpCenterLonF     = 180.
    end if
    ;res@mpMaxLatF        = latn(na) 
    ;res@mpMinLatF        = lats(na)
    ;res@mpMaxLonF        = lonr(na)
    ;res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing 
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;---------set the attribute of the contour-----------------
    res_line = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    if(global) then 
    res_line@gsnAddCyclic     = True
    else
    res_line@gsnAddCyclic     = False
    end if
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = False ;True ;
    res_line@cnLineLabelFontHeightF = 0.01
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.5
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    res_line@gsnContourZeroLineThicknessF = 8
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    reg_line = res_line
    reg_line@cnLineThicknessF = 3.5
    reg_line@cnLineColor      = "green4";"black"
    reg_line@gsnContourZeroLineThicknessF = 0
    reg_line@gsnContourNegLineDashPattern = 5
    reg_line@gsnContourPosLineDashPattern = 0
    reg_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    reg_line@cnLevelSpacingF      := 1     ; set contour interval

    lnres1                      = True
    lnres1@gsLineColor          = "black"
    lnres1@gsLineDashPattern    = 0
    lnres1@gsLineThicknessF     = 5.0

;----------Set res sig-----------
    res_sig                  = True       ; use plot options
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@gsnFrame         = False
    res_sig@gsnDraw          = False
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnLineThicknessF = 2
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "gray50";"darkslategray"
    res_sig@cnFillDotSizeF   = 0.1
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "
    if(global) then 
    res_sig@gsnAddCyclic     = True
    else
    res_sig@gsnAddCyclic     = False
    end if

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait" ;
    if(bar_share) then
    res@lbLabelBarOn      = False;True   ; the label bar
    resp@gsnPanelLabelBar = True
    else
    res@lbLabelBarOn      = True   ; the label bar
    resp@gsnPanelLabelBar = False
    end if
    resp@lbOrientation    = "Vertical" ;"Horizontal" ;
    resp@pmLabelBarWidthF  = 0.03
    resp@pmLabelBarHeightF = 0.4 ;0.27
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.01 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(var_option.eq.-10) then
    minlvl := (/-3.5/)
    maxlvl := (/ 3.5/)
    spacig := (/ 0.5/)
end if
    
if(var_option.eq.-3) then 
    minlvl := (/-7,-0.7/);-8,-8,-8,-8,
    maxlvl := (/ 7, 0.7/); 8, 8, 8, 8,
    spacig := (/ 1, 0.1/); 1, 1, 1, 1,
end if
    
if(var_option.eq.-2) then
    minlvl := (/-0.7,-7,-0.7,-7/);-8,-8,-8,-8,
    maxlvl := (/ 0.7, 7, 0.7, 7/); 8, 8, 8, 8,
    spacig := (/ 0.1, 1, 0.1, 1/); 1, 1, 1, 1,
end if
    
if(var_option.eq.-1) then
    minlvl := (/-0.7,-0.7,-7 /);-8,-8,-8,-8,
    maxlvl := (/ 0.7, 0.7, 7 /); 8, 8, 8, 8,
    spacig := (/ 0.1, 0.1, 1 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.0) then 
    minlvl := (/-2.1,-2.1,-2.1/);
    maxlvl := (/ 2.1, 2.1, 2.1/);
    spacig := (/ 0.3, 0.3, 0.3/);
end if

if(var_option.ge.1) then
if(nnv.eq.0) then   ;dzdt
    minlvl := (/-3.5,-3.5,-3.5,-3.5/);
    maxlvl := (/ 3.5, 3.5, 3.5, 3.5/);
    spacig := (/ 0.5, 0.5, 0.5, 0.5/);
end if
if(nnv.eq.1) then ;dudt
if(var_option.eq.3) then 
    minlvl := (/-21,-21,-21/);-8,-8,-8,-8,
    maxlvl := (/ 21, 21, 21/); 8, 8, 8, 8,
    spacig := (/ 3 , 3 , 3 /); 1, 1, 1, 1,
else
    minlvl := (/-3.5,-3.5,-3.5,-3.5/);-8,-8,-8,-8,
    maxlvl := (/ 3.5, 3.5, 3.5, 3.5/); 8, 8, 8, 8,
    spacig := (/ 0.5, 0.5, 0.5, 0.5/); 1, 1, 1, 1,
end if
end if
if(nnv.eq.2) then ;dTdt
    minlvl := (/-16,-16,-16,-16/);
    maxlvl := (/ 16, 16, 16, 16/);
    spacig := (/  2,  2,  2,  2/);
end if
if(nnv.eq.3) then ;dbdt
    minlvl := (/-24,-24,-24,-24/);
    maxlvl := (/ 24, 24, 24, 24/);
    spacig := (/  3,  3,  3,  3/);
end if
end if

lin_spc = (/1,1,1.5/)

f2 = addfile(path + "/data/ERA-Interim_DJF_uwnd-year.nc","r")
f3 = addfile(filename,"r")
sig   = new((/6,neval/),logical)
eval  = new((/6,neval/),float)
pcvar = new((/6,neval/),float)
do na = 0, narea-1,1
f  = addfile(filename0(na),"r")
ts = f->eof_ts
nyear = dimsizes(ts&year)
eof0 := f->eof0
eof1 := f->eof1
eof2 := f->eof2
eof3 := f->eof3
eof4 := f->eof4
eof5 := f->eof5
    sig(0,:)  = eofunc_north(eof0@pcvar,nyear,False)
    sig(1,:)  = eofunc_north(eof1@pcvar,nyear,False)
    sig(2,:)  = eofunc_north(eof2@pcvar,nyear,False)
    sig(3,:)  = eofunc_north(eof3@pcvar,nyear,False)
    sig(4,:)  = eofunc_north(eof4@pcvar,nyear,False)
    sig(5,:)  = eofunc_north(eof5@pcvar,nyear,False)
    eval(0,:) = eof0@eval
    eval(1,:) = eof1@eval
    eval(2,:) = eof2@eval
    eval(3,:) = eof3@eval
    eval(4,:) = eof4@eval
    eval(5,:) = eof5@eval
    pcvar(0,:) = eof0@pcvar
    pcvar(1,:) = eof1@pcvar
    pcvar(2,:) = eof2@pcvar
    pcvar(3,:) = eof3@pcvar
    pcvar(4,:) = eof4@pcvar
    pcvar(5,:) = eof5@pcvar
ts = ts/conform(ts,sqrt(eval),(/0,1/))
subtitle_color = where(sig,"red","black") 

u := f2->season_ave(:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)})
nlat = dimsizes(u&lat)
nlon = dimsizes(u&lon)
vars := new((/nvar,nyear,nlat,nlon/),double)

if(var_option.le.0.and.var_option.ge.-10)
    do nv = 0, nvar-1,1
    vars(nv,:,:,:) = f3->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    end do
else
    vars  = f3->$var_name(nnv)$(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
if(var_option.eq.-3) then
fbv = addfile(path + "data/" + case + "_DJF_forc_in-year.nc","r")
vars(0,0,:,:) = dim_sum_n_Wrap(vars(:,0,:,:),0) + dim_avg_n_Wrap(fbv->b2(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
vars(1,0,:,:) = dim_avg_n_Wrap(f3->t3x(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) \
               +dim_avg_n_Wrap(f3->t3y(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) 
end if

if(var_option.eq.3) then
fbv = addfile(path + "data/" + case + "_DJF_dzdt_in-year.nc","r")
vars(0,0,:,:) = dim_sum_n_Wrap(vars(0:1,0,:,:),0) + dim_avg_n_Wrap(fbv->$var_name(nnv)$(1,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
vars(1,0,:,:) = dim_avg_n_Wrap(fbv->$var_name(nnv)$(2,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) 
;vars(2,0,:,:) = dim_sum_n_Wrap(vars(0:1,0,:,:),0)
end if
;var_m(1,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)

var  := regCoef_n(ts,vars,2,1)
ureg := regCoef_n(ts,u   ,2,0)
tval := onedtond(var@tval ,(/6,neval,nvar,nlat,nlon/))
df   := onedtond(var@nptxy,(/6,neval,nvar,nlat,nlon/))
utval := onedtond(ureg@tval ,(/6,neval,nlat,nlon/))
udf   := onedtond(ureg@nptxy,(/6,neval,nlat,nlon/))
u    := dim_avg_n_Wrap(u   ,0)
copy_VarCoords(vars(0,0,:,:), var(0,0,0,:,:))
copy_VarCoords(u, ureg(0,0,:,:))
u = mask(u, u.lt.20 , False) 
printVarSummary(var)

b := tval 
b = 0.5
prob  := betainc( df/(df+tval^2), df/2.0, b )
prob   = mask(prob,prob.lt.siglvl,True)
lon_m := conform(prob,vars&lon,4)
prob   = mask(prob,lon_m.eq.0,False)
;var   = mask(var,prob.lt.siglvl,True)
copy_VarCoords(vars(0,0,:,:), prob(0,0,0,:,:))
b := utval 
b = 0.5
uprob  := betainc( udf/(udf+utval^2), udf/2.0, b )
;ureg   = mask(ureg,uprob.lt.siglvl,True)
print(max(tval))
print(min(tval))
print(max(df))
print(min(df))

if(var_option.eq.-6) then 
    var = - var*3600*24
end if

if(var_option.eq.-3) then 
    var(:,:,0,:,:) = var(:,:,0,:,:)*10^11
    var(:,:,1,:,:) = var(:,:,1,:,:)*24*3600
end if

if(var_option.eq.-2) then 
    var(:,:,3,:,:) = var(:,:,3,:,:)*10^11
    var(:,:,2,:,:) = var(:,:,2,:,:)*24*3600
    var(:,:,1,:,:) = var(:,:,1,:,:)*10^11
    var(:,:,0,:,:) = var(:,:,0,:,:)*24*3600
end if

if(var_option.eq.-1) then 
    var(:,:,2,:,:) = var(:,:,2,:,:)*10^11
    var(:,:,0:1,:,:) = var(:,:,0:1,:,:)*24*3600
end if

if(var_option.eq.0) then 
    var = var*24*3600 ;convert units of Q from K/s to K/Day 
end if

if(var_option.ge.1) then
if(nnv.eq.0) then
    var = var*10^3
end if
if(nnv.eq.1) then
    var = var*10^6
end if
if(nnv.eq.2) then
    var = var*10^6
end if
if(nnv.eq.3) then
    var = var*10^12
end if
end if

res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nk = 0, 0, 1
do ne = 0, neval-2, 1
do nv = 0, nvar-1,1
    nev = ne + 1
    np  = nv ;1 var ;2 area
    ;np  = ne*nvar + nv ;1 var ;2 area
    ;np = nk*(neval-1) + ne ;1 eof ;2 kind; 3 area
    if(.not.(global)) then 
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    end if
    res@gsnLeftString   = pre_case(np) + ") " + area(na) + " EOF" + nev ;area(na) + " " + lonl(na) + "-" + lonr(na) + "E"
    res@gsnRightStringFontColor   = subtitle_color(nk,ne)
    if(var_option.ge.1) then 
    res@gsnRightString  = draw_var(nv) + " " + var_name(nnv) + " (" + unit(nnv) + ")  " + eval(nk,ne) + "  " + pcvar(nk,ne)+"%" 
    else
    res@gsnRightString  = draw_var(nv) + " (" + unit(nv) + ")  " + eval(nk,ne) + "  " + pcvar(nk,ne)+"%" 
    end if
    if(nv.gt.10) then
    res@cnLevelSelectionMode  = "AutomaticLevels" ; set manual contour levels
    else
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    end if
    
    plot(np) = gsn_csm_contour_map_ce(wks, var(nk,ne,nv,:,:),res)
    sigl(np) = gsn_csm_contour(wks,prob(nk,ne,nv,:,:),res_sig)
    overlay(plot(np),sigl(np))
   
    reg_line@cnLevelSpacingF      := lin_spc(na)     ; set contour interval
    regl(np) = gsn_csm_contour(wks,ureg(nk,ne,:,:),reg_line)
    overlay(plot(np),regl(np))
   
    line(np) = gsn_csm_contour(wks,u,res_line)
    overlay(plot(np),line(np))
    printMinMax(var(nk,ne,nv,:,:),0)
        
    line1(np) = gsn_add_polyline( wks, plot(np), lonx(:,na), latx(:,na), lnres1 )
end do 
    resp@txString = lev(nl) + "hPa " + area(na); + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
    gsn_panel(wks,plot,(/nvar,1/),resp)
end do 
    ;resp@txString = lev(nl) + "hPa " + draw_var(nv) + " " + area(na); + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
;    resp@txString = var_name(nnv) + " " + area(na); + " (" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)"
;    resp@txString = var_name(nnv) + " " + nk + " EOF" + nev ;+ " " + eval(na,nk,ne) + "  " + pcvar(na,nk,ne)+"%" 
;    gsn_panel(wks,plot,(/3,nvar/),resp)
end do ;EOF kinds
end do ;area 
end do ;var_option
end

