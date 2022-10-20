;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
year_option = True ;calc yearly DJF 
var_option = 1
jet_option  = True ;False;20191205, control the draw of jet stream
region_line = True ;False; draw the line for ME, EA, NA
vect_option = False;True ;
cont_option = False;True ;
bar_share = True ;False;
nl = 4 ;
na = 0 ; "NH","ME","EA","NA"
DJF_option = True ;True ;
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l","(m","(n","(o","(p","(q","(r","(s","(t","(u"/)
;fig_out    = "eps";"x11""eps";
fig_out    = "pdf";"x11""eps";

lev    = (/850,500,300,250,200/) 
if(DJF_option) then
months  = (/"DJF_"/)
months2 = (/"DJF"/)
else
months  = (/""/)
;months2 = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months2 = (/"DJF","MAM","JJA","SON"/)
end if
ntime = dimsizes(months2) 
case   = "ERA-Interim"
;case   = "NCEP1"

do nnv = 1,1,1
    if(var_option.eq.-18) then 
        draw_var = (/"uwnd"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"season_ave"/)  ;use to read variables
        unit     = (/"m/s"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-17) then 
        draw_var = (/"hgt"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"season_ave"/)  ;use to read variables
        unit     = (/"m"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-16) then 
        draw_var = (/"EGR"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"season_ave"/)  ;use to read variables
        unit     = (/"1/day"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-15) then 
        draw_var = (/"T"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"season_ave"/)  ;use to read variables
        unit     = (/"~S~o~N~C"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-14) then 
        draw_var := (/"Ah","Al"/)  ;,"fva","adv","A","cur" this is the variable used to draw, can used as center string
        var_name := (/"Ah","Al"/)  ;,"fva","adv","A","cur" use to read variables
        unit     := "10~S~-5~N~" + (/"m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-13) then 
        draw_var := (/"-fua","adv","A","cur"/)  ; this is the variable used to draw, can used as center string
        var_name := (/"fva","adv","A","cur"/)  ; use to read variables
        unit     := "10~S~-5~N~" + (/"m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-12) then 
        draw_var = (/"vort"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"vort"/)  ;use to read variables
        unit     = (/"10~S~-6~N~ s~S~-1~N~"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-11) then 
        draw_var = (/"dTdy"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"dTdy"/)  ;use to read variables
        unit     = (/"10~S~-6~N~K/m"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-10) then 
        var_calc = (/"intAPE","intEKE","variV","variZ"/)
        draw_var = (/"all"   ,"hgh"   ,"low"   /)  ;this is the variable used to draw, can used as center string
        var_name = (/"intEKE","hghEKE","lowEKE"/)  ;use to read variables
        unit     = (/"10~S~4~N~J","10~S~4~N~J","m~S~2~N~s~S~-2~N~","10~S~4~N~m~S~4~N~s~S~-4~N~"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-6) then 
        draw_var := (/"T3x","T3y","Tadvect"/)  ;Tadvect
        var_name := (/"t3x","t3y","b3y"/)  ;this is the variable used to draw, can used as center string
        unit     := (/"K/day","k/day","k/day"/)  ;this can be used as right string 
        select := (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.-5) then 
        draw_var := (/"vort","b2","total"/)  ;,"b3y" this is the variable used to draw, can used as center string
        ;draw_var := (/"vort_x","vort_y","vort"/)  ;,"b3y" this is the variable used to draw, can used as center string
        var_name := (/"b1x","b1y","b2"/)  ;,"b3y" use to read variables
        unit     := "10~S~-11~N~" + (/"s~S~-2~N~","s~S~-2~N~","s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
        select := (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.-4) then 
        draw_var := (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
        var_name := (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
        unit     := "10~S~-11~N~" + (/"s~S~-2~N~","s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
    end if
    if(var_option.eq.-3) then 
        draw_var := (/"Qte_hgh","Fte_hgh"/)  ;this is the variable used to draw, can used as center string
        var_name := (/"Qeddh","Ah"/)  ;use to read variables
        unit     := (/"K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
        select := (/0,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.-2) then 
        draw_var := (/"Qte_low","Fte_low"/)  ;this is the variable used to draw, can used as center string
        var_name := (/"Qeddl","Al"/)  ;use to read variables
        unit     := (/"K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
        select := (/1,3/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.-1) then 
        draw_var := (/"Qd"   ,"Qte"  ,"Fte"/)  ;this is the variable used to draw, can used as center string
        var_name := (/"Qd"   ,"Qeddy","A"  /)  ;use to read variables
        unit     := (/"K/Day","K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
        select := (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.0) then 
        draw_var = (/"Qd"   ,"Qte","total"/)  ;this is the variable used to draw, can used as center string
        var_name = (/"Qd"   ,"Qeddy"/)  ;use to read variables
        unit     = (/"K/Day","K/Day","K/day"/)  ;this can be used as right string 
    end if
    if(var_option.ge.1.and.var_option.le.6) then
    var_name := (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
    unit     := (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
    ;nnv = 1 ;used to read dzdt induced by different forcing
    if(var_option.eq.1) then 
        draw_var := (/"Qd","Qte","Fte"/)  ;this is the variable used to draw, can used as center string
        select := (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.2) then 
        draw_var := (/"Qte_hgh","Fte_hgh","Qte_low","Fte_low"/)  ;this is the variable used to draw, can used as center string
        select := (/0,2,1,3/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.3) then 
        draw_var := (/"Qte_hgh","Fte_hgh"/)  ;this is the variable used to draw, can used as center string
        select := (/0,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.4) then 
        draw_var := (/"b1","b2","Tadvect"/)  ;,"total" this is the variable used to draw, can used as center string
        select := (/0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
    if(var_option.eq.5) then 
        ;draw_var := (/"Tx","Ty","Tadvect"/)  ;,"b3y" this is the variable used to draw, can used as center string
        ;select := (/0,1,2/)
        ;draw_var := (/"vort_x","vort_y","vort"/)  ;,"b3y" this is the variable used to draw, can used as center string
        draw_var := (/"vort","b2","total"/)  ;,"b3y" this is the variable used to draw, can used as center string
        select := (/0,1,2/);(/2,3,1/) 
    end if
    if(var_option.eq.6) then 
        draw_var := (/"Qd","Qte","Fte","b1","b2","b3x","b3y","turn","b1z","b4","total","dudt"/)  ;this is the variable used to draw, can used as center string
        select := (/0,1,2,0,1,0,1,0,1,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
    end if
end if
if(var_option.eq.13) then
    var_name = (/"fva","adv","A","Ah","Al","cur"/)  ; use to read variables
    unit     = "10~S~-5~N~" + (/"m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~"/)  ;this can be used as right string 
;    nnv = 3 ;used to read dzdt induced by different forcing
    if(nnv.eq.0) then 
        draw_var := (/"fv","-fvg","fva"/)
        select   := (/0,1/)
    end if
    if(nnv.eq.1) then 
        draw_var := (/"-ududx","-vdudy","-wdudz","adv"/)  ;this is the variable used to draw, can used as center string
        select   := (/0,1,2/)
    end if
    if(nnv.ge.2.and.nnv.le.4) then 
        draw_var := (/"-duudx","-duvdy","-duwdz","A"/)  ;,"uv/a","-uw/a" this is the variable used to draw, can used as center string
        select   := (/0,1,2/)
    end if
    if(nnv.eq.5) then 
        draw_var := (/"uv","uw","cur"/)  ;this is the variable used to draw, can used as center string
        select   := (/0,1/)
    end if
end if
if(var_option.eq.14) then
    var_name = (/"fva","adv","A","Ah","Al","cur"/)  ; use to read variables
    unit     = "10~S~-5~N~" + (/"m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~","m*s~S~-2~N~"/)  ;this can be used as right string 
;    nnv = 3 ;used to read dzdt induced by different forcing
    if(nnv.eq.0) then 
        draw_var := (/"-fu","fug","-fua"/)
        select   := (/0,1/)
    end if
    if(nnv.eq.1) then 
        draw_var := (/"-udvdx","-vdvdy","-wdvdz","adv"/)  ;this is the variable used to draw, can used as center string
        select   := (/0,1,2/)
    end if
    if(nnv.ge.2.and.nnv.le.4) then 
        draw_var := (/"-duvdx","-dvvdy","-dvwdz","-uu/a","-vw/a","A"/)  ; this is the variable used to draw, can used as center string
        select   := (/0,1,2,3,4/)
    end if
    if(nnv.eq.5) then 
        draw_var := (/"uu","vw","cur"/)  ;this is the variable used to draw, can used as center string
        select   := (/0,1/)
    end if
end if
if(var_option.eq.20) then 
    draw_var = (/"E_all" ,"E_low" ,"E_hgh"  /)  ;this is the variable used to draw, can used as center string
    var_name = (/"ex_all","ey_all","ez_all" /)  ;
    unit     = (/"m2/s2" ,"m2/s2" ,"10~S~-2~N~Pa*m/s~S~2~N~"/)  ;this can be used as right string 
    divg_option = False ;True
end if
if(var_option.eq.21) then 
    draw_var = (/"W_all" ,"W_low" ,"W_hgh"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"wafx"  ,"wafy"  ,"wafz" /)  ;
    unit     = (/"m2/s2" ,"m2/s2" ,"10~S~-5~N~m/s~S~2~N~"/)  ;this can be used as right string 
    divg_option = True ;False;
end if
if(var_option.eq.22) then 
    draw_var = (/"steady"/)  ;this is the variable used to draw, can used as center string
    var_name = (/"wafx"  ,"wafy"  ,"edd_sf" /)  ;
    unit     = (/"m2/s2" ,"m2/s2" ,"m~S~2~N~/s"/)  ;this can be used as right string 
    divg_option = False;True ;
end if

;area: east Asia
area = (/"NH","ME","EA","NA","TP"/)
lats = (/ 10 , 0  ,  0 ,  0 , 10 /);10 ;-20
latn = (/ 60 , 89 , 89 , 89 , 50 /) ;55 ;60
lonl = (/ 0  , 0  ,105 , 250, 40 /) ;30 ;20
lonr = (/360 , 70 ,175 , 320,130 /) ;150;220
lonx = (/ lonl, lonr, lonr, lonl, lonl /)
latx = (/ lats, lats, latn, latn, lats /)
lat_spacing = 20
lon_spacing = 30
density = 1.1 ;km/m3
g = 9.8 ;m2/s2

;the path of the file
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
;gtopfile = "/home/ys17-19/data-observation/NCEP_NCAR/hgt.sfc.nc"
if(var_option.le.-16)
    filename = path + "data/" + case + "_" + months + draw_var(0) + "-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + draw_var(0)
end if
if(var_option.eq.-15)
    filename = path + "data/" + case + "_" + months + "air-year.nc"
    contname = path + "data/" + case + "_"+months+"hgt-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.eq.-13.or.var_option.eq.-14)
    filename = path + "data/" + case + "_" + months + "meridional_momentum-year.nc"
    ;filename = path + "data/" + case + "_" + months + "zonal_momentum-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.eq.-10)
    filename = path + "data/" + case + "_" + months + var_calc(nnv) + "-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_calc(nnv)
end if
if(var_option.eq.-5.or.var_option.eq.-6)
    filename = path + "data/" + case + "_"+months+"forc_in3-year.nc"
    contname = path + "data/" + case + "_"+months+"dzdt_in3-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + lev(nl) + var_name(0)
end if
if(var_option.eq.-4.or.var_option.eq.-12)
    filename = path + "data/" + case + "_"+months+"forc_in-year.nc"
    contname = path + "data/" + case + "_"+months+"dzdt_in-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.le.-2.and.var_option.ge.-3)
    filename = path + "data/" + case + "_"+months+"forc_te-year.nc"
    contname = path + "data/" + case + "_"+months+"dzdt_te-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.le.0.and.var_option.ge.-1)
    filename = path + "data/" + case + "_"+months+"forc_ex-year.nc"
    contname = path + "data/" + case + "_"+months+"dzdt_ex-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.eq.1)
    filename = path + "data/" + case + "_"+months+"dzdt_ex-year.nc"
    fig_name = path + "fig/1p2x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.ge.2.and.var_option.le.3)
    filename = path + "data/" + case + "_"+months+"dzdt_te-year.nc"
    fig_name = path + "fig/1p2x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.4)
    filename = path + "data/" + case + "_"+months+"dzdt_in-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.5)
    filename = path + "data/" + case + "_DJF_dzdt_b1xb1yb2-year.nc"
    ;filename = path + "data/" + case + "_DJF_dzdt_in3-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.6)
    filename  = path + "data/" + case + "_DJF_dzdt_ex-year.nc"
    filename1 = path + "data/" + case + "_DJF_dzdt_in-year.nc"
    filename2 = path + "data/" + case + "_DJF_dzdt_in3-year.nc"
    filename3 = path + "data/" + case + "_DJF_dzdt_re-year.nc"
    filename4 = path + "data/" + case + "_DJF_forc_re-year.nc"
    fig_name  = path + "fig/1p3x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.13)
    filename = path + "data/" + case + "_" + months + "zonal_momentum_flux-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.14)
    filename = path + "data/" + case + "_" + months + "meridional_momentum_flux-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(nnv) + "_" + draw_var(0)
end if
if(var_option.eq.20)
    filename = path + "data/" + case + "_" + months + "e_vector-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.eq.21)
    filename = path + "data/" + case + "_" + months + "TN01-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if
if(var_option.eq.22)
    filename = path + "data/" + case + "_" + months + "steady_TN01-year.nc"
    fig_name = path + "fig/1p3x3_horizontal_" + var_name(0)
end if

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)})) ;hgt for NCEP1,unit is gpm ;z for ERA-I, unit is m2/s2
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nlat   = dimsizes(phis&latitude )  ;
nlon   = dimsizes(phis&longitude)  ;
nvar   = dimsizes(draw_var)
var_m  := new((/nvar,ntime,nlat,nlon/),double)
if(var_option.ge.20) then
var_m  := new((/3,nvar,ntime,nlat,nlon/),double)
end if
if(cont_option) then
var_c  := new((/nvar,ntime,nlat,nlon/),double)
end if

;the number of the panel plot
    plot     = new((/12/),graphic)
    line     = new((/12/),graphic)
    vect     = new((/12/),graphic)
    line1    = new((/12/),graphic)
    line2    = new((/12/),graphic)
    line3    = new((/12/),graphic)
    gtop     = new((/12/),graphic)
    cont     = new((/12/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.01
   vcref_font     = 0.005
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.5
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    if(var_option.eq.-10.or.var_option.eq.-15.or.var_option.le.-16) then 
    gsn_define_colormap(wks,"precip_11lev")
    ;gsn_define_colormap(wks,"precip3_16lev")
    ;gsn_define_colormap(wks,"precip2_17lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    else
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    end if
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red

;do nl = 0,4, 1
;================================================================
;read data
;====================================================================
if(jet_option) then
    f = addfile(path + "data/" + case + "_" + months + "uwnd-year.nc","r")
    u_m = new((/ntime,nlat,nlon/),float)
    if(DJF_option) then
    u_m(0,:,:) = dim_avg_n_Wrap(f->season_ave(:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    else
    u_m = dim_avg_n_Wrap(f->month_ave(:,:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    end if
    u_m = mask(u_m, u_m.lt.25 , False) 
    printVarSummary(u_m)
end if

if(vect_option) then
    u  = new((/ntime,nlat,nlon/),float)
    v  = new((/ntime,nlat,nlon/),float)
    f1 = addfile(path + "data/" + case + "_" + months + "uwnd-year.nc","r")
    f2 = addfile(path + "data/" + case + "_" + months + "vwnd-year.nc","r")
    if(DJF_option) then
    u(0,:,:) = dim_avg_n_Wrap(f1->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    v(0,:,:) = dim_avg_n_Wrap(f2->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    else
    u = dim_avg_n_Wrap(f1->season_ave(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    v = dim_avg_n_Wrap(f2->season_ave(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    end if
end if

if(var_option.eq.-11) then
    f = addfile(path + "data/" + case + "_" + months + "air-year.nc","r")
    t = new((/ntime,nlat,nlon/),float)
    if(DJF_option) then
    t(0,:,:) = dim_avg_n_Wrap(f->season_ave(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    else
    t = dim_avg_n_Wrap(f->season_ave(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    end if
    
    opt    = 0     ;used by center_finite_diff_n, no meanging 
    loncyc = True ;used by center_finite_diff_n
    cyclic = False ;used by center_finite_diff_n
    ;means NOT treated as cyclic in differenced coordinates
    ;the end values will use a one-sided difference scheme for the end points. 
    a  = 6378388 ;the radius of earth, m
    pi = atan(1.0)*4
    lat = t&lat
    lat = lat*pi/180.0
    coslat = conform(var_m(0,:,:,:),cos(lat),1)
;    var_m(0,:,:,:) = tofloat(center_finite_diff_n(t*coslat,lat,cyclic,opt,1)/a/coslat)
    var_m(0,:,:,:) = tofloat(center_finite_diff_n(t,lat,cyclic,opt,1)/a)
    var_m = var_m*1000000
    copy_VarMeta(t,var_m(0,:,:,:))
end if

if((var_option.le.0.and.var_option.ge.-10).or.var_option.le.-12)
f = addfile(filename,"r")
if(var_option.eq.-10.and.nnv.le.1) then 
    do nv = 0, nvar-1,1
    if(DJF_option) then
    var_m(nv,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m(nv,:,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,:,{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    end do
else
    do nv = 0, nvar-1,1
    if(DJF_option) then
    var_m(nv,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m(nv,:,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    end do
    var_m(nvar-1,:,:,:) = dim_sum_n_Wrap(var_m(0:nvar-2,:,:,:),0)
    if(cont_option) then 
    fc = addfile(contname,"r")
    if(var_option.ge.-6) then 
        if (nnv.eq.0) then 
            var_c(:,0,:,:) = dim_avg_n_Wrap(fc->dzdt(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
            var_c = var_c*10^3
        else
            var_c(:,0,:,:) = dim_avg_n_Wrap(fc->dudt(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
            var_c = var_c*10^6
        end if
        var_c(nvar-1,:,:,:) = dim_sum_n_Wrap(var_c(0:nvar-2,:,:,:),0)
    else
        if(DJF_option) then
        var_c(0,0,:,:) = dim_avg_n_Wrap(fc->$var_name(0)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
        else
        var_c(0,:,:,:) = dim_avg_n_Wrap(fc->$var_name(0)$(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
        end if
        var_c = var_c/9.8
    end if 
    end if 
end if
end if

if(var_option.ge.1.and.var_option.le.5)
f = addfile(filename,"r")
    if(DJF_option) then
    var_m(:,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m  = dim_avg_n_Wrap(f->$var_name(nnv)$(select,:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end if
if(var_option.eq.6)
    f = addfile(filename,"r")
    var_m(0:2,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(select(0:2),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    f = addfile(filename1,"r")
    var_m(3:4,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(select(3:4),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    f = addfile(filename2,"r")
    var_m(5:6,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(select(5:6),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    f = addfile(filename3,"r")
    var_m(7:9,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(select(7:9),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    var_m(nvar-2,:,:,:) = dim_sum_n_Wrap(var_m(0:nvar-3,:,:,:),0)
    f = addfile(filename4,"r")
    var_m(nvar-1,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    copy_VarMeta(phis,var_m(0,0,:,:))
end if
if(var_option.ge.13.and.var_option.lt.20)
f = addfile(filename,"r")
    if(DJF_option) then
    var_m(0:nvar-2,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(:,select,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) 
    else
    var_m(0:nvar-2,:,:,:)  = dim_avg_n_Wrap(f->$var_name(nnv)$(:,select,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) 
    end if
    var_m(nvar-1,:,:,:) = dim_sum_n_Wrap(var_m(0:nvar-2,:,:,:),0)
end if

if(var_option.eq.-5) then
;fbv = addfile(path + "data/" + case + "_" + months + "forc_in-year.nc","r")
;var_m(1,0,:,:) = var_m(1,0,:,:) + dim_avg_n_Wrap(f->b2(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
;var_m(2,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)
var_m(0,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)
var_m(1,0,:,:) = dim_avg_n_Wrap(f->b2(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
var_m(2,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)
if(cont_option) then 
    fbv = addfile(path + "data/" + case + "_" + months + "dzdt_in-year.nc","r")
    var_c(1,0,:,:) = var_c(1,0,:,:) + dim_avg_n_Wrap(fbv->dudt(1,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)*10^6
    var_c(2,:,:,:) = dim_sum_n_Wrap(var_c(0:1,:,:,:),0)
end if
end if

if(var_option.eq.5) then
;fbv = addfile(path + "data/" + case + "_" + months + "dzdt_in-year.nc","r")
;var_m(1,0,:,:) = var_m(1,0,:,:) + dim_avg_n_Wrap(fbv->$var_name(nnv)$(1,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
;var_m(1,:,:,:) = dim_sum_n_Wrap(var_m(1:2,:,:,:),0)
var_m(0,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)
var_m(1,:,:,:) = var_m(2,:,:,:)
var_m(2,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)
end if
;var_m(1,:,:,:) = dim_sum_n_Wrap(var_m(0:1,:,:,:),0)

if(var_option.ge.20) then 
f = addfile(filename,"r")
do nv = 0,2,1
    if(DJF_option) then
    if(var_option.eq.22) then 
    var_m(nv,0,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    else
    var_m(nv,:,0,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    else
    var_m(nv,:,:,:,:)  = dim_avg_n_Wrap(f->$var_name(nv)$(:,:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    end if
end do
if(divg_option) then 
    var_m(2,:,:,:,:) = uv2dv_cfd( var_m(0,:,:,:,:), var_m(1,:,:,:,:),var_m&lat, var_m&lon, 3 )*100000
    minlvl := (/-3.5,-3.5,-3.5/); -1.4,-1.4,-1.4-3.5,-3.5,-3.5-7,-7,-7
    maxlvl := (/ 3.5, 3.5, 3.5/);  1.4, 1.4, 1.4 3.5, 3.5, 3.5 7, 7, 7
    spacig := (/ 0.5, 0.5, 0.5/);  0.2, 0.2, 0.2 0.5, 0.5, 0.5 1, 1, 1
else
    if(var_option.eq.22) then ;(1.0/1000000.0)*
    var_m(2,0,0,:,:) = dim_avg_n_Wrap(f->$var_name(2)$(:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;get an array (ncase,nmonth,nlat,nlon)
    var_m(2,:,0,:,:) = -100.0*dim_avg_n_Wrap(f->$var_name(2)$(:,:,{925},{lats(na):latn(na)},{lonl(na):lonr(na)}),1) ;get an array (ncase,nmonth,nlat,nlon)
    end if
    minlvl := (/-21,-21,-21/);-3.5,-3.5,-3.54.2,-4.2,-4.2
    maxlvl := (/ 21, 21, 21/); 3.5, 3.5, 3.54.2, 4.2, 4.2
    spacig := (/ 3 , 3 , 3 /); 0.5, 0.5, 0.50.6, 0.6, 0.6
end if
vcref  := (/50,150,150/) ;(/925,850,700,500,200/)
;vcref  := (/150,150,150/) ;(/925,850,700,500,200/)
;vcref  := (/50,50,50/) ;(/925,850,700,500,200/)
end if
printVarSummary(var_m)
printMinMax(var_m,True)

if(var_option.eq.-17) then 
    var_m = var_m/9.8
end if

if(var_option.eq.-16) then 
    var_m = var_m*24*3600
end if

if(var_option.eq.-15) then 
    var_m = var_m-273.15
end if

if(var_option.le.-13.and.var_option.ge.-14) then 
    var_m = var_m*100000
end if

if(var_option.eq.-12) then
    pi = atan(1.0)*4
    lat  = (var_m&lat)*pi/180.0
;    var_m = var_m + conform(var_m,2*(2*pi/24.0/3600.0)*sin(lat),2)
    var_m = var_m*1000000
end if

if(var_option.eq.-10) then 
if(nnv.le.1) then 
    var_m = var_m/10000
    minlvl := (/10 ,10 ,10 /);-8,-8,-8,-8,
    maxlvl := (/130,130,130/); 8, 8, 8, 8,
    spacig := (/8  ,8  ,8  /); 1, 1, 1, 1,
end if
if(nnv.eq.2) then 
    minlvl := (/20 ,10 ,10 /);-8,-8,-8,-8,
    maxlvl := (/320,160,160/); 8, 8, 8, 8,
    spacig := (/30 ,15 ,15 /); 1, 1, 1, 1,
end if
if(nnv.eq.3) then 
    var_m = var_m/10000
    minlvl := (/10 ,0  ,10 /);-8,-8,-8,-8,
    maxlvl := (/310,105,310/); 8, 8, 8, 8,
    spacig := (/20 ,7 ,20  /); 1, 1, 1, 1,
end if
end if

if(var_option.eq.-6) then 
    var_m = - var_m*3600*24
end if

if(var_option.eq.-4.or.var_option.eq.-5) then 
    var_m = var_m*10^11
end if

if(var_option.le.-2.and.var_option.ge.-3)
    var_m(1,:,:,:) = var_m(1,:,:,:)*10^11
    var_m(0,:,:,:) = var_m(0,:,:,:)*24*3600
end if

if(var_option.eq.-1) then 
    var_m(2,:,:,:)   = var_m(2,:,:,:)*10^11
    var_m(0:1,:,:,:) = var_m(0:1,:,:,:)*24*3600
end if

if(var_option.eq.0) then 
    var_m = var_m*24*3600 ;convert units of Q from K/s to K/Day 
end if

if(var_option.ge.1.and.var_option.le.6) then
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

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
;    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.6
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 150.
    res@mpMaxLatF        = latn(na) 
    res@mpMinLatF        = lats(na)
    if(na.eq.0) then
    res@gsnAddCyclic     = True;False
    res@mpMaxLonF        = 330 ;lonr(na)
    res@mpMinLonF        = -30 ;lonl(na)
    else
    res@gsnAddCyclic     = False;True;
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    end if
    res@gsnMajorLatSpacing = lat_spacing 
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
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
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
;    res@cnLineThicknessF = 1.0

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
    
    ;res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.2
    res@pmLabelBarParallelPosF       = 0.5
    ;res@pmLabelBarWidthF  = 0.7
    ;res@pmLabelBarHeightF = 0.1

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    
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

;--------set the attribute of vector--------------
    resvc                  = True
    resvc@gsnDraw          = False
    resvc@gsnFrame         = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMonoLineArrowColor       = True
    resvc@vcLineArrowColor           = "gray20"
    resvc@vcMinDistanceF             = 0.015
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcLineArrowHeadMinSizeF   = 0.005
    resvc@vcLineArrowHeadMaxSizeF   = 0.005
    ;resvc@vcRefMagnitudeF            = 40.0             ; define vector ref mag
    ;resvc@vcRefAnnoString2           = "40m/s"
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
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

;---------set the attribute of the contour-----------------
    res_line = True
    res_line@gsnDraw          = False
    res_line@gsnFrame         = False
    res_line@gsnAddCyclic     = True;False
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True ;False
    res_line@cnLineLabelFontHeightF = 0.005
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "blueviolet";"green4";"black"
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 2
    res_line@gsnContourPosLineDashPattern = 0
    
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@lbLabelBarOn     = False   ; the label bar
    res_line@gsnCenterString  = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    res_cont = True
    res_cont@gsnDraw          = False
    res_cont@gsnFrame         = False
    res_cont@gsnAddCyclic     = True;False
    res_cont@cnFillOn         = False
    res_cont@cnLinesOn        = True
    res_cont@cnLineLabelsOn   = False ;True ;
    res_cont@cnLineLabelFontHeightF = 0.005
    res_cont@cnInfoLabelOn    = False
    res_cont@cnLineThicknessF = 3.5
    if(jet_option) then 
    res_cont@cnLineColor      = "green4";"black";"blueviolet";
    else
    res_cont@cnLineColor      = "blueviolet";"green4";"black";
    end if
    res_cont@gsnContourZeroLineThicknessF = 0
    res_cont@gsnContourNegLineDashPattern = 7
    res_cont@gsnContourPosLineDashPattern = 0
  
    levels = new((/10/),integer)
    levels(0:4) = (/-300,-250,-150,-80,-20/)
    levels(5:9) = -levels(0:4:-1)
    ;res_cont@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    ;res_cont@cnLevels             = levels       ; set contour interval
    res_cont@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_cont@cnLevelSpacingF      = 1        ; set contour interval
    res_cont@lbLabelBarOn     = False   ; the label bar
    res_cont@gsnCenterString  = ""
    res_cont@gsnLeftString    = " "
    res_cont@gsnRightString   = " "

    if(na.eq.0) then 
    res_cont@gsnAddCyclic     = True;False
    res_line@gsnAddCyclic     = True;False
    else
    res_cont@gsnAddCyclic     = False
    res_line@gsnAddCyclic     = False
    end if

;--------set the attribute of polygon--------------
    lnres1                      = True
    lnres1@gsLineColor          = "black"
    lnres1@gsLineDashPattern    = 0
    lnres1@gsLineThicknessF     = 5.0
    lnres2       = lnres1
    lnres2@gsLineDashPattern    = 2
    lnres3       = lnres1
    lnres3@gsLineDashPattern    = 16

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    if(bar_share) then
    res@lbLabelBarOn      = False;True   ; the label bar
    resp@gsnPanelLabelBar = True
    else
    res@lbLabelBarOn      = True   ; the label bar
    resp@gsnPanelLabelBar = False
    end if
    resp@lbOrientation     = "Horizontal" ;"Vertical" ;
    resp@pmLabelBarWidthF  = 0.8          ;0.05       ;
    resp@pmLabelBarHeightF = 0.06         ;0.8        ;
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(var_option.le.-6) then 
    vcref = (/10,20,25,30,40/) ;(/925,850,700,500,200/)
end if

if(var_option.eq.-18) then ;Eady growth rate
    minlvl := (/20/);
    maxlvl := (/70/);
    spacig := (/5 /);
end if

if(var_option.eq.-17) then ;Eady growth rate
    minlvl := (/2820/);
    maxlvl := (/3140/);
    spacig := (/20  /);
end if

if(var_option.eq.-16) then ;Eady growth rate
    minlvl := (/ 0 /);
    maxlvl := (/1.6/);
    spacig := (/0.1/);
end if

if(var_option.eq.-15) then ;T ;(/925,850,700,500,250/)
if(nl.eq.0) then ;850 
    minlvl := (/-16/);
    maxlvl := (/ 16/);
    spacig := (/ 2 /);
end if
if(nl.eq.1) then ;500 
    minlvl := (/-26/);
    maxlvl := (/-10/);
    spacig := (/ 1 /);
end if
if(nl.eq.2) then ;300
    minlvl := (/-51/);
    maxlvl := (/-35/);
    spacig := (/ 1 /);
end if
if(nl.eq.3) then  ;250
    minlvl := (/-58/);
    maxlvl := (/-42/);
    spacig := (/ 1 /);
end if
if(nl.eq.4) then ;200
    minlvl := (/-60/);-64
    maxlvl := (/-52/);-48
    spacig := (/ 1 /); 1
end if
end if

if(var_option.eq.-14) then
    minlvl := (/-7,-7/);
    maxlvl := (/ 7, 7/);
    spacig := (/ 1, 1/);
end if
if(var_option.eq.-13) then
    minlvl := (/-35,-35,-14,-14/);-8,-8,-8,-8,
    maxlvl := (/ 35, 35, 14, 14/); 8, 8, 8, 8,
    spacig := (/ 5 , 5 , 2 , 2 /); 1, 1, 1, 1,
end if
if(var_option.eq.-12) then
    minlvl := (/-49 /)
    maxlvl := (/ 49 /)
    spacig := (/ 7  /)
end if
if(var_option.eq.-11) then
    minlvl := (/-14 /)
    maxlvl := (/ 14 /)
    spacig := (/ 2 /)
end if
if(var_option.eq.-6) then
    minlvl := (/-10.5,-10.5,-10.5/);-8,-8,-8,-8,
    maxlvl := (/ 10.5, 10.5, 10.5/); 8, 8, 8, 8,
    spacig := (/ 1.5, 1.5, 1.5/); 1.5, 1.5, 1.5, 1.5,
end if
if(var_option.eq.-4.or.var_option.eq.-5) then 
    minlvl := (/-70,-70,-70,-70/);-14,-14,-14,-14
    maxlvl := (/ 70, 70, 70, 70/); 14, 14, 14, 14
    spacig := (/ 10 , 10 , 10 , 10 /); 2 , 2 , 2 , 2 
end if
if(var_option.le.-2.and.var_option.ge.-3)
    minlvl := (/-2.1,-14/);-8,-8,-8,-8,
    maxlvl := (/ 2.1, 14/); 8, 8, 8, 8,
    spacig := (/ 0.3, 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.-1) then
    minlvl := (/-2.1,-2.1,-14/);-8,-8,-8,-8,
    maxlvl := (/ 2.1, 2.1, 14/); 8, 8, 8, 8,
    spacig := (/ 0.3, 0.3, 2 /); 1, 1, 1, 1,
end if
    
if(var_option.eq.0) then 
    minlvl := (/-2.1,-2.1,-2.1/);
    maxlvl := (/ 2.1, 2.1, 2.1/);
    spacig := (/ 0.3, 0.3, 0.3/);
end if

if(var_option.ge.1.and.var_option.le.6) then
    minlvl = (/-3.5,-70,-16,-24/);-35,35
    maxlvl = (/ 3.5, 70, 16, 24/); 35,35
    spacig = (/ 0.5, 10,  2,  3/);  5,5 
end if

if(var_option.ge.13.and.var_option.lt.20) then ;dbdt
    minlvl = (/-70,-35,-14,-7,-7,-7/);-8,-8,-8,-8,
    maxlvl = (/ 70, 35, 14, 7, 7, 7/); 8, 8, 8, 8,
    spacig = (/ 10, 5 , 2 , 1, 1, 1/); 1, 1, 1, 1,
end if

    ;res@cnFillPalette = "BlueDarkRed18"
    ;blue and red, has 18 colors,2-10 is blue
    
    if(.not.(var_option.eq.-10.or.var_option.eq.-15.or.var_option.le.-16)) then 
    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
    end if
do nt = 0, ntime-1,1
do nv = 0,nvar-1,1
;do ni = 0, 3, 1
;do np = 0, 2, 1 
    ;nt = ni*3 + np
    np = nv
    if(var_option.ge.20) then 
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa " + draw_var(nv) + " WAF (" + unit(0) + ")"
    if(divg_option) then 
    res@gsnRightString  = "divg (" + unit(2) + ")"
    else
    res@gsnRightString  = var_name(2) + " (" + unit(2) + ")"
    end if
    end if
    if(var_option.gt.0.and.var_option.lt.20) then 
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa " + draw_var(nv)
    res@gsnRightString  = months2(nt) + " " + var_name(nnv) + " (" + unit(nnv) + ")"
    end if
    if(var_option.le.0) then
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa" 
    res@gsnRightString  = months2(nt) + " " + draw_var(nv) + " (" + unit(nv) + ")"
    if(var_option.eq.-11) then 
    res@gsnLeftString   = pre_case(np) + ") " + months2(nt) 
    end if
    end if
    if(var_option.eq.-10) then 
    res@gsnLeftString   = pre_case(np) + ") " + lev(nl) + "hPa " + var_calc(nnv)
    res@gsnRightString  = months2(nt) + " " + draw_var(nv) + " (" + unit(nnv) + ")"
    end if
    
    ;res@gsnSpreadColors       = False              ; turn of the full color map
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    if(var_option.ge.1)
    res@cnLevelSpacingF       := spacig(nnv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nnv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nnv)       ; set min contour level
    else
    res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    end if
    if(var_option.ge.20)
    plot(np) = gsn_csm_contour_map_ce(wks,var_m(2,nv,nt,:,:),res)
    resvc@vcRefMagnitudeF      = vcref(nv)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nv)+"m~S~2~N~/s~S~2~N~"
    vect(np) = gsn_csm_vector(wks,var_m(0,nv,nt,:,:),var_m(1,nv,nt,:,:),resvc)
    overlay(plot(np),vect(np))
    printMinMax(var_m(0,nv,nt,:,:),0)
    printMinMax(var_m(1,nv,nt,:,:),0)
    printMinMax(var_m(2,nv,nt,:,:),0)
    else
    plot(np) = gsn_csm_contour_map_ce(wks,var_m(nv,nt,:,:),res)
    printMinMax(var_m(nv,nt,:,:),0)
    end if

    gtop(np) = gsn_csm_contour(wks,phis,res_gtop)
    overlay(plot(np),gtop(np))
    
    if(vect_option) then 
    resvc@vcRefMagnitudeF      = vcref(nl)             ; define vector ref mag
    resvc@vcRefAnnoString2     = vcref(nl)+"m/s"
    vect(np) = gsn_csm_vector(wks,u(nt,:,:),v(nt,:,:),resvc)
    overlay(plot(np),vect(np))
    end if

    if(cont_option) then 
    cont(np) = gsn_csm_contour(wks,var_c(nv,nt,:,:),res_cont)
    overlay(plot(np),cont(np))
    printMinMax(var_c(nv,nt,:,:),0)
    end if
    
    if(region_line) then 
        line1(np) = gsn_add_polyline( wks, plot(np), lonx(:,1), latx(:,1), lnres1 )
        line2(np) = gsn_add_polyline( wks, plot(np), lonx(:,2), latx(:,2), lnres2 )
        line3(np) = gsn_add_polyline( wks, plot(np), lonx(:,3), latx(:,3), lnres3 )
    end if
    
    if(jet_option) then 
    line(np) = gsn_csm_contour(wks,u_m(nt,:,:),res_line)
    overlay(plot(np),line(np))
    end if
end do 
;    gsn_panel(wks,plot,(/3,1/),resp)
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/3,1/),resp)
;     gsn_panel(wks,plot,(/nvar,1/),resp)
end do ;lev 
;     gsn_panel(wks,plot,(/2,2/),resp)
;end do 
end

