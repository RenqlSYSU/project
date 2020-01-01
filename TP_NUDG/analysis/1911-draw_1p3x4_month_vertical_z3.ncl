;*********************************************************************
;Draw the shaded map of monthly or seasonal average precipitation in east asia. 
;And there are three cases,1AMIP-CTRL,2CTRL-NDUG24h,3CTRL-NUDG6h 
;then do the significance testing, the dot mean the difference is significance
;                                                     by Ql Ren
;                                                    2019/02/01
;*********************************************************************************
begin
path    = "/home/ys17-19/renql/project/TP_NUDG"
path1   = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h","NUDG24h"/)
    case2   = (/"CTRL","NG6h","NG24h"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
else
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
end if

draw_var = (/"Z3" ,"dzdt" ,"Ahgt","U","T","dTdy","dudt"/)  ;this is the variable used to draw, can used as center string
unit     = (/"gpm","10~S~-4~N~m~S~2~N~/s~S~3~N~","gpm","m/s","~S~o~N~C","10~S~-6~N~m/s","10~S~-5~N~m/s~S~2~N~"/)    ;this can be used as right string 
draw_varname = (/"hgt","dzdt","Ahgt","U","T","-dTdy","dudt"/)  ;this is the variable used to draw, can used as center string
lev   = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels ,hPa U
lev_u = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels ,hPa U
nv = 6
diff = True;False
jet_option = True
jetfile  = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_U.nc"

var_option = 2
if(var_option.eq.0) then 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
imonth = ispan(0,3,1)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
end if
if(var_option.eq.1) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
imonth = ispan(0,11,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
end if
if(var_option.eq.2) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"JF","MA","MJ","JA"/)
;months   = (/"MJ","JA"/)
imonth = ispan(0,3,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_MJ_JA_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
end if
if(var_option.eq.3) then 
var_name = draw_var(nv) ;this is the var_name in data file
months   = (/"JF","MA","MJ","JA"/)
imonth = ispan(0,3,1);(/4,5,6,7/)
month_s = 1
month_e = 8
filein  = path1+(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_MJ_JA_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
lev  := (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
lev_u:= (/1000,925,850,700,600,500,400,300,200,100/)
years    = ispan(1979,2005,1)   ;the calculate time
end if

;pre_case   = (/"(m","(n","(o","(p","(e","(f","(g","(h","(i","(j","(k","(l"/)
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    = "eps";"pdf";"x11"

;area: east Asia
area = (/"TP","EA","NH"/)
lats = (/15  ,15  , 0  /);20 ;-20
latn = (/55  ,55  , 89 /) ;55 ;60
lonl = (/70  ,110 , 10 /) ;30 ;20
lonr = (/105 ,130 ,350 /) ;150;220
lat_spacing = 15
lon_spacing = 30
na = 1

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case1)
nmonth = dimsizes(imonth)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nlev   = dimsizes(lev)
nlev2  = dimsizes(lev_u)
var    = new((/ncase,nmonth,nlev ,nlat,nlon/),float)
var_u  = new((/ncase,nmonth,nlev2,nlat/),float)

;the number of the panel plot
    plot   = new((/ncase*nmonth/),graphic)
    line   = new((/ncase*nmonth/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 1.8
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
printVarSummary(var)
if(var_option.eq.3) then 
    f    = addfile(filein(0)+years(0)+".daily."+var_name+".nc","r")
    time = cd_calendar(f->time,0)
    select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
    ntime = dimsizes(select_time)
    nyear = dimsizes(years)
    var0 = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
    
    do nf = 0 , ncase-1,1
    do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name+".nc","r")
    var0(nf,nt,:,:,:,:) = f->$var_name$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    end do 
    end do
    var0 := dim_avg_n_Wrap(var0,1)
    
    do nm = 0, nmonth-1, 1
    mons = 2*nm + 1
    mone = 2*nm + 2
    select_time1 := ind(time(:,1).ge.mons.and.time(:,1).le.mone)
    var(:,nm,:,:,:) = dim_avg_n_Wrap(var0(:,select_time1,:,:,:),1)
    end do
else
    do nf = 0 ,ncase-1, 1
    print(filename(nf) + draw_var(nv) + ".nc")
    f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
    var(nf,:,:,:,:)  = f->$var_name(0)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end do
    if(jet_option) then 
    do nf = 0 ,ncase-1, 1
    f = addfile(jetfile(nf),"r")
    if(var_option.eq.0) then
    u := f->season_ave(:,{lev_u},{lats(na):latn(na)},{lonl(na):lonr(na)})
    else
    u := f->month_ave(:,{lev_u},{lats(na):latn(na)},{lonl(na):lonr(na)})
    end if
    var_u(nf,:,:,:) = dim_avg_n_Wrap(u,3)
    end do
    end if
end if
printVarSummary(var)
var := dim_avg_n_Wrap(var,4)

if(jet_option) then 
    var_u = mask(var_u,var_u.ge.0,True)
   ; copy_VarMeta(vars(0,:,:,0),var_u(0,0,:,:))
    printVarSummary(var_u)
end if

if(nv.eq.1) then
    var = var*9.8*10000 ;dzdt,transform unit from m/s to m2/s3
end if

;if(nv.eq.3) then
;    u = mask(var, var.lt.0 , False) 
;    copy_VarMeta(var,u)
;end if

if(nv.eq.4) then
    var = var-273.15 ;T,transform unit from K to C
end if

if(nv.eq.5) then
    var = -var*10^6 ;transform dTdy to -dTdy
end if

if(nv.eq.6) then
    var = var*10^5 ;dudt
end if

if(diff) then 
do nf = 0 ,ncase-2,1
case2   = (/"CTRL","NG6h-CTRL","NG24h-CTRL"/)
var(nf+1,:,:,:)= var(nf+1,:,:,:) - var(0,:,:,:)   ;get the difference between CTRL and NUDG
end do
end if
printVarSummary(var)

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
    
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-32,-24,-16,-10,-5,-1,-0.5,-0.1, 0, 0.1, 0.5, 1, 5, 10, 16, 24, 32/)

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
    
    res@lbLabelBarOn         = True; the label bar False ;
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0
    res@pmLabelBarParallelPosF       = 0.5

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
    res_line = res
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.0
    res_line@cnLineColor      = "black"
   
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = False;True
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.04
    resp@lbLabelFontHeightF = 0.012
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(nv.eq.0) then 
    minlvl = (/  100,-35,-35/);
    maxlvl = (/16600, 35, 35/);
    spacig = (/ 1100, 5 , 5 /);
end if

if(nv.eq.1) then 
    minlvl = (/-4 ,-0.8,-0.8/);
    maxlvl = (/ 4 , 0.8, 0.8/);
    spacig = (/0.5, 0.1, 0.1/);
end if

if(nv.eq.2) then 
    minlvl = (/-80,-24,-24/);
    maxlvl = (/ 80, 24, 24/);
    spacig = (/ 10, 3 , 3 /);
end if

if(nv.eq.3) then 
    minlvl = (/-28,-7,-7/);
    maxlvl = (/ 28, 7, 7/);
    spacig = (/ 4 , 1, 1/);
end if

if(nv.eq.4) then 
    minlvl = (/-40,-2.1,-2.1/);
    maxlvl = (/ 20, 2.1, 2.1/);
    spacig = (/ 4 , 0.3, 0.3/);
end if

if(nv.eq.5) then 
    minlvl = (/ 0 ,-2.4,-2.4/);
    maxlvl = (/ 60, 2.4, 2.4/);
    spacig = (/ 4 , 0.3, 0.3/);
end if

if(nv.eq.6) then 
if(diff) then
    minlvl = (/-21,-7,-7/);
    maxlvl = (/ 21, 7, 7/);
    spacig = (/ 3 , 1, 1/);
else
    minlvl = (/-21,-21,-21/);
    maxlvl = (/ 21, 21, 21/);
    spacig = (/ 3 , 3 , 3 /);
end if
end if

    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
do nc = 0, ncase-1,1
do nm = 0,nmonth-1,1
;if((nc.gt.0).or.(nv.eq.1).or.(nv.eq.2).or.(nv.eq.3)) then
;    res@cnFillPalette = "BlueDarkRed18"
;   ;blue and red, has 18 colors,2-10 is blue
;else
;   res@cnFillPalette = "precip3_16lev"
;   ;blue and yellow and red, has 17 colors,2-18,2 is white
;end if 
        np = nmonth*nc + nm
        ;np = nm
    ;    np = ncase*nm + nc
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(imonth(nm)) 
        res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_varname(nv) ;+ " (" + unit(nv) + ")"
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        plot(np) = gsn_csm_pres_hgt(wks,var(nc,nm,:,:),res)
        if(jet_option) then 
        line(np) = gsn_csm_pres_hgt(wks,var_u(nc,nm,:,:),res_line)
        overlay(plot(np),line(np))
        end if
        print(max(var(nc,nm,:,:)))
        print(min(var(nc,nm,:,:)))
    end do 
    ;gsn_panel(wks,plot,(/1,4/),resp)
end do ;for the month 

    ;gsn_panel(wks,plot,(/nmonth,ncase/),resp)
    gsn_panel(wks,plot,(/ncase,nmonth/),resp)
;end do
end

