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
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h"/) ;,"NUDG24h"
    case2   = (/"CTRL","NG6h"/) ;,"NG24h"
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
    case1   = (/"F2000","NGTR"/)       ;NGTR    
    case2   = (/"F2000","NGTR"/) ;NGTR
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 30
end if

draw_var = (/"Z3","dzdt","Ahgt","U","T", "dTdy","dudt","4q1","V","EKE"/)  ;this is the variable used to draw, can used as center string
unit     = (/"gpm","10~S~-4~N~m~S~2~N~/s~S~3~N~","gpm","m/s","~S~o~N~C","10~S~-6~N~m/s","10~S~-5~N~m/s~S~2~N~","K/day","m/s","J"/)    ;this can be used as right string 
draw_varname = (/"hgt","dzdt","Ahgt","U","T","-dTdy","dudt","4q1","V","EKE"/)  ;this is the variable used to draw, can used as center string
;lev   = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels ,hPa U
lev   = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;13 levels ,hPa , 4q1
;lev_u = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels ,hPa U
lev_u = lev
nv = 9
diff = False;True ;
perc = True
jet_option = True ;False 
barshare = True ;False;
barvert  = False;True
siglvl   = 0.05

var_option = 1
if(var_option.eq.0) then 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
imonth = ispan(0,3,1)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
jetfile    = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
end if
if(var_option.eq.1) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Ann","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
imonth = ispan(0,11,1);(/4,5,6,7/)
annual = True
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
jetfile    = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
end if
if(var_option.eq.2) then 
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"JF","MA","MJ","JA"/)
;months   = (/"MJ","JA"/)
imonth = ispan(2,3,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-vertical_section_MJ_JA_" + draw_var(nv) + "_" + var_name(0) + "-3c" 
jetfile    = path + "/analysis/mdata/" + case1 + "-Clim_2month_ave_U.nc"
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
;lev  := (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;hPa
years    = ispan(1979,2005,1)   ;the calculate time
end if
if(var_option.eq.4) then 
var_name = (/"EKE","EKE_vari"/) ;this is the var_name in data file
months   = (/"ANN"/)
imonth = ispan(0,0,1);(/4,5,6,7/)
filename = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_" ;+ draw_var(nv) + ".nc"
fig_name = path + "/analysis/fig/" + case1(0) + "-level_lon_" + draw_var(nv) + "_" + var_name(0)
jetfile    = path + "/analysis/mdata/" + case1 + "-Clim_month_ave_U.nc"
end if

pre_case   = (/"(g","(h","(o","(p","(e","(f","(g","(h","(i","(j","(k","(l"/)
;pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    = "pdf";"x11"
;fig_out    = "eps";"pdf";"x11"

;area: east Asia
area = (/"TP","EA","NH","CA"/)
lats = (/ 35 , 15 , 0  , 15 /);20 ;-20
latn = (/ 35 , 55 , 89 , 55 /) ;55 ;60
lonl = (/  0 ,110 , 1  ,0   /) ;30 ;20
lonr = (/180 ,130 ,359 ,70  /) ;150;220
lon_spacing = 30
na = 0

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase  = dimsizes(case1)
nmonth = dimsizes(imonth)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nlev   = dimsizes(lev)
nlev2  = dimsizes(lev_u)
if(nv.eq.7) then 
var    = new((/ncase,4,nmonth,nlev ,nlat,nlon/),float)
else
var    = new((/ncase,nmonth,nlev ,nlat,nlon/),float)
vari   = new((/ncase,nmonth,nlev ,nlat,nlon/),float)
;ps     = new((/ncase,nmonth,nlat,nlon/),float)
end if
var_u  = new((/ncase,nmonth,nlev2,nlon/),float)
phis1  := 1013.25*(1-dim_avg_n_Wrap(phis,0)*0.0065/288.15)^5.25145
phis1   = where(phis1.le.1000, phis1, 1000)
ypoly = new((/nlon+4/),float)
xpoly = new((/nlon+4/),float)
ypoly(0:nlon-1)    = tofloat(phis1)
ypoly(nlon:nlon+3) = (/tofloat(phis1(nlon-1)),1000,1000,tofloat(phis1(0))/)
xpoly(0:nlon-1)    = tofloat(phis&lon)
xpoly(nlon:nlon+3) = (/tofloat(lonr(na)),tofloat(lonr(na)),tofloat(lonl(na)),tofloat(lonl(na))/)
delete(phis1)

;the number of the panel plot
    plot   = new((/ncase*nmonth/),graphic)
    line   = new((/ncase*nmonth/),graphic)
    core   = new((/ncase*nmonth/),graphic)
    sig    = new((/ncase*nmonth/),graphic)
    topo   = new((/ncase*nmonth/),graphic)

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
    var := dim_avg_n_Wrap(var,4)
end if

if(var_option.le.2) then 
    if(nv.eq.7) then 
        do nf = 0 ,ncase-1, 1
        print(filename(nf) + draw_var(nv) + ".nc")
        f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
        var(nf,:,:,:,:,:)  = f->$var_name(0)$(:,imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end do
        var := dim_avg_n_Wrap(var,4)
    else
        do nf = 0 ,ncase-1, 1
        print(filename(nf) + draw_var(nv) + ".nc")
        f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
        var(nf,:,:,:,:)  = f->$var_name(0)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        vari(nf,:,:,:,:) = f->$var_name(1)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        if(annual) then 
        var(nf,0,:,:,:)  = f->EKE(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        vari(nf,0,:,:,:) = f->EKE_vari(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
        end if
        end do
        var := dim_avg_n_Wrap(var,3)
        vari := dim_avg_n_Wrap(vari,3)
    end if
end if

if(var_option.eq.4) then 
    do nf = 0 ,ncase-1, 1
    print(filename(nf) + draw_var(nv) + ".nc")
    f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
    var(nf,0,:,:,:) = f->$var_name(0)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    vari(nf,0,:,:,:) = f->$var_name(1)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end do
    var  := dim_avg_n_Wrap(var,3)
    vari := dim_avg_n_Wrap(vari,3)
end if
    
;vari  = where(abs(vari).le.0.0000001, 0.01, vari)
prob  = ttest(var(0,:,:,:),vari(0,:,:,:),nsample,var(1,:,:,:),vari(1,:,:,:),nsample,True,False)
prob  = mask(prob, prob.lt.siglvl, True)
copy_VarMeta(var(0,:,:,:),prob)
printVarSummary(var)

if(jet_option) then 
    do nf = 0 ,ncase-1, 1
    f = addfile(jetfile(nf),"r")
    if(var_option.eq.0) then
    u := f->season_ave(imonth,{lev_u},{lats(na):latn(na)},{lonl(na):lonr(na)})
    else
    u := f->month_ave(imonth,{lev_u},{lats(na):latn(na)},{lonl(na):lonr(na)})
    end if
    var_u(nf,:,:,:) = dim_avg_n_Wrap(u,2)
    end do
    
    if(annual) then 
    var_u(:,0,:,:) = dim_avg_n_Wrap(var_u,1)
    end if
    
    copy_VarMeta(u(0,:,0,:),var_u(0,0,:,:))
    printVarSummary(var_u)
    corez = new((/ncase,nmonth/),float)
    corex = new((/ncase,nmonth/),float)
    do nc = 0, ncase-1, 1
    do nt = 0, nmonth-1, 1
        term = ndtooned(var_u(nc,nt,:,:))
        indices = ind_resolve(maxind(term),(/nlev2,nlon/))
        corez(nc,nt) = tofloat(var_u&lev(indices(0,0)))
        corex(nc,nt) = tofloat(var_u&lon(indices(0,1)))
    end do 
    end do 
    printVarSummary(indices)
    print(corez)
    print(corex)
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

if(nv.eq.7) then
    var = var*3600*24 ;4q1,k/s to k/day
end if

if(diff) then 
case2(1)   = case1(1) + "-" + case1(0) ;,"NG24h-CTRL"
do nf = 0 ,ncase-2,1
if(nv.eq.7) then 
var(nf+1,:,:,:,:)= var(nf+1,:,:,:,:) - var(0,:,:,:,:)   ;get the difference between CTRL and NUDG
else
if(perc) then 
var(nf+1,:,:,:)= var(nf+1,:,:,:)*100/var(0,:,:,:)-100   ;get the difference between CTRL and NUDG
else
var(nf+1,:,:,:)= var(nf+1,:,:,:) - var(0,:,:,:)   ;get the difference between CTRL and NUDG
end if
end if
end do
end if
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)

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
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
;    res@cnLineThicknessF = 1.0
    
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-32,-24,-16,-10,-5,-1,-0.5,-0.1, 0, 0.1, 0.5, 1, 5, 10, 16, 24, 32/)

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
    res@tmXBTickSpacingF   = 30
    res@tmXBTickStartF     = lonl(na)
    res@tmXBTickEndF       = lonr(na)

;----------Set res sig-----------
    res_sig                  = res       ; use plot options
    res_sig@cnFillOn         = True       ; Fill contours
    res_sig@cnLinesOn        = False;True      ; turn off contour lines
    res_sig@gsnFrame         = False
    res_sig@gsnDraw          = False
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels
    res_sig@cnLineThicknessF = 2
    ;res_sig@cnLineColor     = "black"
    res_sig@cnInfoLabelOn    = False
    res_sig@cnMonoFillPattern= True;False
    res_sig@cnFillPattern    = 11
    res_sig@cnMonoFillScale  = True
    res_sig@cnFillScaleF     = 0.5
    res_sig@cnMonoFillColor  = True
    res_sig@cnFillColor      = "darkslategray"
    res_sig@cnFillDotSizeF   = 0.25
    res_sig@lbLabelBarOn     = False;True   ; turn off the label bar
    ;res_sig@cnLevelSelectionMode ="ExplicitLevels"
    ;res_sig@cnLevels             = (/3000,1000/)
    res_sig@gsnLeftString    = " "
    res_sig@gsnRightString   = " "

;---------set the attribute of the contour-----------------
    res_line = res
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "black"
    res_line@gsnContourZeroLineThicknessF = 8
    res_line@gsnContourNegLineDashPattern = 2
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

;used to draw the jet core
    polyres                   = True
    polyres@gsMarkerIndex     = 5        ; polymarker style
    polyres@gsMarkerSizeF     = 8.          ; polymarker size 
    polyres@gsMarkerColor     = "black"
    polyres@gsMarkerThicknessF  = 5.0 

;used to draw the black terrain
    pres                  = True
    pres@gsFillColor      = "gray"
    
    res@cnMissingValFillColor = "gray"
    
    resp  = True    
    resp@gsnMaximize   = False;True
    resp@gsnPaperOrientation = "landscaoe";"portrait"
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
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
    res@pmLabelBarHeightF = 0.1
    res@pmLabelBarOrthogonalPosF     = 0.1 ;y direction
    res@pmLabelBarParallelPosF       = 0.5
    end if
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFont        = font_index              ; make labels smaller
    resp@lbLabelFontHeightF = 0.015
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

if(nv.eq.7) then 
    minlvl = (/-2.1,-2.1,-2.1/);-0.7,-0.7
    maxlvl = (/ 2.1, 2.1, 2.1/); 0.7, 0.7
    spacig = (/ 0.3, 0.3, 0.3/); 0.1, 0.1
end if

if(nv.eq.8) then 
    minlvl = (/-4 ,-0.8,-0.8/);
    maxlvl = (/ 4 , 0.8, 0.8/);
    spacig = (/0.5, 0.1, 0.1/);
end if

if(nv.eq.9) then 
if(diff) then
    minlvl = (/ 10,-70, 10/);
    maxlvl = (/160, 70,160/);
    spacig = (/ 10, 10, 10/);
else
    minlvl = (/ 10, 10/);
    maxlvl = (/160,160/);
    spacig = (/ 10, 10/);
end if
end if

if(nv.eq.7) then
draw_varname2 = (/"QRL","QRS","DTCOND","DTV","total"/) 
res@cnFillPalette = "BlueDarkRed18"
res@cnFillColors  := (/0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17/)
;blue and red, has 18 colors,2-10 is blue
do nv2 = 0,4,1
do nc = 0, ncase-1,1
do nm = 0,nmonth-1,1
        np = nm
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(imonth(nm)) 
        res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_varname2(nv2) + " (" + unit(nv) + ")"
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        if(nv2.eq.4) then 
        plot(np) = gsn_csm_pres_hgt(wks,dim_sum_n_Wrap(var(nc,:,nm,:,:),0),res)
        else
        plot(np) = gsn_csm_pres_hgt(wks,var(nc,nv2,nm,:,:),res)
        print(max(var(nc,nv2,nm,:,:)))
        print(min(var(nc,nv2,nm,:,:)))
        end if
        if(jet_option) then 
        line(np) = gsn_csm_pres_hgt(wks,var_u(nc,nm,:,:),res_line)
        overlay(plot(np),line(np))
        end if
    end do 
    gsn_panel(wks,plot,(/3,4/),resp)
end do ;for the month 
end do ;for the month 

else
do nc = 0, ncase-1,1
do nm = 0,0,1 ;nmonth-1,1
if((nc.le.0).or.(nv.eq.9)) then
    res@cnFillPalette = "precip3_16lev"
    res@cnFillColors  := (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
    ;blue and yellow and red, has 17 colors,2-18,2 is white
else
    res@cnFillPalette = "BlueDarkRed18"
    res@cnFillColors  := (/0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17/)
    ;blue and red, has 18 colors,2-10 is blue
end if 
        ;np = nmonth*nc + nm
         np = nc
    ;    np = ncase*nm + nc
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(imonth(nm)) 
        res@gsnRightString  = lats(na) + "N " + draw_varname(nv) + " (" + unit(nv) + ")";+ "-" + latn(na) 
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        plot(np) = gsn_csm_pres_hgt(wks,var(nc,nm,:,:),res)
        
        if(nc.gt.0.and.diff) then 
        sig(np) = gsn_csm_pres_hgt(wks,prob(nm,:,:),res_sig)
        overlay(plot(np),sig(np))
        end if
        
        if(jet_option) then 
        line(np) = gsn_csm_pres_hgt(wks,var_u(nc,nm,:,:),res_line)
        overlay(plot(np),line(np))
        ;core(np) = gsn_add_polymarker(wks,plot(np),corey(1,nm),corez(1,nm),polyres)
        end if
        topo(np) = gsn_add_polygon(wks, plot(np), xpoly, ypoly, pres)
        print(max(var(nc,nm,:,:)))
        print(min(var(nc,nm,:,:)))
    end do 
;    gsn_panel(wks,plot,(/3,4/),resp)
end do ;for the month 
    
;    gsn_panel(wks,plot,(/nmonth,ncase/),resp)
    ;gsn_panel(wks,plot,(/ncase,nmonth/),resp)
    gsn_panel(wks,plot,(/1,2/),resp)
end if
end

